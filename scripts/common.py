"""Small shared I/O and reference checks used by the fetcher and analyzer."""
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote

MISSING = {'', '.', './.', '.|.', 'NA', 'N/A', 'NAN', 'NONE', 'NULL', '?', 'N'}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def empty_output(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise ValueError(f'Output directory must be empty: {path}')
    path.mkdir(parents=True, exist_ok=True)
    return path


def promoter_window(start, end, strand, length):
    if strand not in ['+', '-'] or start < 1 or end < start or length < 1:
        raise ValueError('Invalid gene interval, strand, or promoter length')
    if strand == '+':
        if start == 1:
            raise ValueError('No upstream chromosome sequence before position 1')
        return max(1, start - length), start - 1
    return end + 1, end + length


def cds_features(path, transcript):
    opener = gzip.open if str(path).endswith('.gz') else open
    segments, locus = [], set()
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            c = line.rstrip().split('\t')
            if len(c) != 9 or c[2] != 'CDS':
                continue
            attrs = dict(item.split('=', 1) for item in c[8].split(';') if '=' in item)
            if transcript not in [unquote(x) for x in attrs.get('Parent', '').split(',')]:
                continue
            if c[7] not in ['0', '1', '2']:
                raise ValueError('Missing CDS phase; cannot validate automatic translation')
            locus.add((c[0], c[6]))
            segments.append((int(c[3]), int(c[4]), int(c[7])))
    if not segments or len(locus) != 1:
        raise ValueError(f'Expected one CDS locus for {transcript}')
    chromosome, strand = locus.pop()
    if strand not in ['+', '-']:
        raise ValueError('Unknown CDS strand')
    segments.sort(reverse=strand == '-')
    positions, phase_ok = [], True
    for start, end, phase in segments:
        phase_ok &= phase == ((3 - len(positions) % 3) % 3)
        positions.extend(range(start, end + 1) if strand == '+' else range(end, start - 1, -1))
    if len(positions) != len(set(positions)):
        raise ValueError('Overlapping CDS features')
    return chromosome, strand, positions, bool(phase_ok)


def reconstruct_proteins(haps, sites, sequence, strand):
    from Bio.Seq import Seq
    sequence = sequence.upper()
    if not sequence or len(sequence) % 3 or set(sequence) - set('ACGT'):
        raise ValueError('Reference CDS must be complete, unambiguous, and divisible by three')
    reference_peptide = str(Seq(sequence).translate())
    if not reference_peptide.startswith('M') or '*' in reference_peptide.rstrip('*'):
        raise ValueError('Reference CDS has no start methionine or has internal stops')
    comp = str.maketrans('ACGT', 'TGCA')
    result = {}
    used = set()
    for row in sites.to_dict('records'):
        idx = int(row['cds_position']) - 1
        if idx < 0 or idx >= len(sequence) or idx in used:
            raise ValueError('Invalid or duplicate CDS position')
        used.add(idx)
        alleles = str(row['alleles']).upper().split('/')
        if len(alleles) < 2 or any(len(a) != 1 or a not in 'ACGT' for a in alleles):
            raise ValueError('Automatic protein reconstruction supports explicit SNP alleles only')
        coding = [a.translate(comp) if strand == '-' else a for a in alleles]
        if sequence[idx] not in coding:
            raise ValueError(f'Reference mismatch at {row["site_id"]}')
    for _, hap in haps.iterrows():
        mutated = list(sequence)
        for row in sites.to_dict('records'):
            allele = str(hap[row['site_id']]).upper()
            if allele not in str(row['alleles']).upper().split('/'):
                raise ValueError('Missing or unsupported haplotype allele; no reference imputation')
            mutated[int(row['cds_position']) - 1] = allele.translate(comp) if strand == '-' else allele
        result[hap['haplotype_id']] = str(Seq(''.join(mutated)).translate())
    return result, reference_peptide
