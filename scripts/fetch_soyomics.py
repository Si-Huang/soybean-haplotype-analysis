"""Download one verified SoyZH13 gene into a portable, auditable bundle."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

import pandas as pd
import requests

from common import cds_features, empty_output, promoter_window, sha256, write_json
from style import POPULATIONS

BASE = 'https://ngdc.cncb.ac.cn/soyomics'
GFF_URL = 'https://download.cncb.ac.cn/public/soyomics/gene/ZH13.gene.gff.gz'


class Archive:
    def __init__(self, directory):
        self.directory = directory
        self.records = []
        self.session = requests.Session()
        self.session.headers['User-Agent'] = 'soybean-haplotype-analysis/1.0'

    def get(self, name, url, data=None, binary=False):
        response = self.session.request('POST' if data is not None else 'GET', url,
                                       data=data, timeout=(15, 180))
        path = self.directory / name
        path.write_bytes(response.content)
        self.records.append({'url': url, 'method': response.request.method,
                             'parameters': data, 'file': name, 'status': response.status_code,
                             'fetched_utc': datetime.now(timezone.utc).isoformat(),
                             'sha256': sha256(path)})
        write_json(self.directory / 'provenance.json', self.records)
        response.raise_for_status()
        return response.content if binary else response.text

    def json(self, name, url, data=None):
        return json.loads(self.get(name, url, data))


def scalar(value):
    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError('Expected scalar or length-one JSON list')
        return value[0]
    return value


def finite_filter(value, lower, upper):
    try:
        return lower <= float(value) <= upper
    except (ValueError, TypeError):
        return False


def sequence_from_json(payload):
    text = payload['data']
    text = ''.join(text) if isinstance(text, list) else text
    text = re.sub(r'<br\s*/?>', '\n', text, flags=re.I)
    seq = ''.join(line.strip() for line in text.splitlines() if not line.startswith('>')).upper()
    if not seq or set(seq) - set('ACGTN'):
        raise ValueError('Sequence endpoint returned unexpected FASTA content')
    return seq


def fetch(args):
    if not re.fullmatch(r'SoyZH13_\d{2}G\d+', args.gene):
        raise ValueError('Resolve the gene identifier to ZH13.v2 first; do not rename ID prefixes')
    if not 0 <= args.maf <= 0.5 or not 0 <= args.max_missing <= 1 or args.min_count < 1:
        raise ValueError('Invalid filtering thresholds')
    out = empty_output(args.out)
    raw = out / 'raw'
    raw.mkdir()
    api = Archive(raw)
    print(f'Fetching {args.gene} ({args.role})', flush=True)
    search = api.json('gene_search.json', f'{BASE}/search/genes?id={args.gene}&draw=1&start=0&length=10')
    genes = [g for g in search['data'] if g['gene_id'] == args.gene]
    if len(genes) != 1:
        raise ValueError('Gene identifier did not resolve uniquely')
    gene = genes[0]
    html = api.get('gene.html', f'{BASE}/genes/{args.gene}')
    strand_match = re.search(r'<strong>Strand:</strong></td>\s*<td[^>]*>([+-])</td>', html)
    if not strand_match:
        raise ValueError('Could not verify strand from gene record')
    strand = strand_match[1]
    start, end = int(gene['gene_start']), int(gene['gene_end'])
    chromosome = gene['chr']
    transcript = args.transcript or args.gene + '.m1'
    manifest = {'gene_id': args.gene, 'symbol': args.symbol or args.gene,
                'role': args.role, 'assembly': 'ZH13.v2', 'strand': strand,
                'gene_start': start, 'gene_end': end, 'maf': args.maf,
                'max_missing': args.max_missing, 'min_haplotype_copies': args.min_count,
                'scope_note': 'Frequencies describe the returned HapSnap subset; phase provenance is service-derived.',
                'created_utc': datetime.now(timezone.utc).isoformat()}
    mapping = {}
    if args.role == 'target':
        lo, hi = promoter_window(start, end, strand, args.promoter_length)
        manifest.update({'promoter_length': args.promoter_length,
                         'boundary_note': 'Upstream of annotated gene boundary; chromosome-end length not independently verified.'})
    else:
        if not transcript.startswith(args.gene + '.'):
            raise ValueError('Transcript must belong to the resolved gene')
        manifest['transcript'] = transcript
        if args.gff:
            gff = Path(args.gff)
            manifest['gff_sha256'] = sha256(gff)
        else:
            gff = raw / 'ZH13.gene.gff.gz'
            api.get(gff.name, GFF_URL, binary=True)
        cds_chr, cds_strand, positions, phase_ok = cds_features(gff, transcript)
        if (cds_chr, cds_strand) != (chromosome, strand):
            raise ValueError('GFF locus does not agree with gene query; check genome version')
        mapping = {p: i + 1 for i, p in enumerate(positions)}
        lo, hi = min(positions), max(positions)
        payload = api.json('reference_cds.json', BASE + '/plumber/seq_fetch',
                           {'strain': 'ZH13', 'gtype': 'trans', 'region_kw': f'{chromosome}:{lo}-{hi}',
                            'rna_type': 'mRNA', 'gene_kw': transcript, 'seq_type': 'cds'})
        sequence = sequence_from_json(payload)
        manifest['cds_phase_valid'] = phase_ok and len(sequence) == len(positions)
        (out / 'reference_cds.fasta').write_text(f'>{transcript}\n{sequence}\n', encoding='utf-8')
    manifest.update({'chromosome': chromosome, 'region_start': lo, 'region_end': hi})
    variants = api.json('variants_all.json', f'{BASE}/genome_variation_list2?chr={chromosome}&start_position={lo}&end_position={hi}')
    if not isinstance(variants, list):
        raise ValueError('Variant endpoint schema changed; inspect archived response')
    selected = [v for v in variants if finite_filter(v.get('maf'), args.maf, 0.5)
                and finite_filter(v.get('miss'), 0, args.max_missing)
                and (args.role == 'target' or int(v['position']) in mapping)]
    selected.sort(key=lambda v: int(v['position']))
    if len({v['position'] for v in selected}) != len(selected):
        raise ValueError('Multiple records per genomic position; normalize locally before haplotyping')
    site_rows = []
    for v in selected:
        p = int(v['position'])
        site_rows.append({'site_id': str(p), 'position': p, 'var_id': v.get('var_id', ''),
                          'relative_position': p - start if strand == '+' else end - p,
                          'cds_position': mapping.get(p, ''), 'alleles': v.get('alleles', ''),
                          'maf': v['maf'], 'missing': v['miss'], 'consequence': v.get('consequence', '')})
    columns = ['site_id', 'position', 'var_id', 'relative_position', 'cds_position', 'alleles', 'maf', 'missing', 'consequence']
    pd.DataFrame(site_rows, columns=columns).to_csv(out / 'sites_requested.csv', index=False)
    manifest['requested_site_count'] = len(site_rows)
    hap_rows, counts, assignments, retained = [], [], [], []
    status = 'no_filtered_variants'
    if selected:
        body = {'varList': json.dumps([f'{chromosome}:{v["position"]}' for v in selected]),
                'imputation': '0', 'individual_miss': '0.1', 'site_miss': str(args.max_missing),
                'maf': str(args.maf), 'heter': '1', 'populations': ','.join(POPULATIONS),
                'count': str(args.min_count)}
        write_json(raw / 'hapsnap_request.json', body)
        hap = api.json('hapsnap_response.json', BASE + '/plumber/ldhap', body)
        if hap.get('error') and ''.join(map(str, hap['error'] if isinstance(hap['error'], list) else [hap['error']])):
            raise ValueError(f'HapSnap error: {hap["error"]}')
        data = hap.get('data', [])
        if not isinstance(data, list):
            raise ValueError('Unexpected HapSnap data schema')
        status = 'ok' if data else 'no_haplotypes_returned'
        if data:
            retained = [r for r in site_rows if r['site_id'] in data[0]]
            if not retained:
                raise ValueError('HapSnap returned types without retained input sites')
            source_map = {}
            for index, row in enumerate(data, 1):
                source_id = str(scalar(row['Hap_ID']))
                if source_id in source_map:
                    raise ValueError('Duplicate HapSnap source ID')
                hid = f'H{index}'
                source_map[source_id] = hid
                hap_rows.append({'haplotype_id': hid, 'source_id': source_id,
                                 **{r['site_id']: scalar(row[r['site_id']]) for r in retained}})
                for pop in POPULATIONS:
                    counts.append({'haplotype_id': hid, 'population': pop, 'count': int(scalar(row[pop]))})
            for row in hap.get('sample_data', []):
                source_id = str(scalar(row['Hap_ID']))
                if source_id not in source_map:
                    raise ValueError('Unknown type in sample_data')
                for pop in POPULATIONS:
                    values = row[pop] if isinstance(row[pop], list) else [row[pop]]
                    for sample_id in values:
                        if str(sample_id) not in ['', '0', 'None']:
                            assignments.append({'sample_id': str(sample_id),
                                                'haplotype_id': source_map[source_id], 'population': pop})
    sites = [r['site_id'] for r in retained]
    pd.DataFrame(retained, columns=columns).to_csv(out / 'sites.csv', index=False)
    pd.DataFrame([{'site_id': r['site_id'], 'retained': r['site_id'] in sites} for r in site_rows],
                 columns=['site_id', 'retained']).to_csv(out / 'site_retention.csv', index=False)
    pd.DataFrame(hap_rows, columns=['haplotype_id', 'source_id'] + sites).to_csv(out / 'haplotypes.csv', index=False)
    pd.DataFrame(counts, columns=['haplotype_id', 'population', 'count']).to_csv(out / 'counts.csv', index=False)
    if assignments:
        pd.DataFrame(assignments).to_csv(out / 'assignments.csv', index=False)
        metadata = []
        sample_ids = sorted({r['sample_id'] for r in assignments})
        for i in range(0, len(sample_ids), 250):
            rows = api.json(f'metadata_{i // 250}.json', BASE + '/haplotype/strain_list?strain_list',
                            {'strain_list': ','.join(sample_ids[i:i + 250])})
            if not isinstance(rows, list):
                raise ValueError('Unexpected sample metadata schema')
            for r in rows:
                metadata.append({'sample_id': r['accession'], 'accession': r.get('strainName', ''),
                                 'country': r.get('country', ''), 'region': r.get('region', ''),
                                 'longitude': r.get('longitude', ''), 'latitude': r.get('latitude', '')})
        pd.DataFrame(metadata).to_csv(out / 'metadata.csv', index=False)
    manifest.update({'site_columns': sites, 'status': status})
    write_json(out / 'manifest.json', manifest)
    print(f'Wrote {len(hap_rows)} types at {len(sites)} retained sites: {out}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gene', required=True)
    p.add_argument('--role', required=True, choices=['target', 'tf'])
    p.add_argument('--symbol')
    p.add_argument('--transcript')
    p.add_argument('--gff')
    p.add_argument('--out', required=True)
    p.add_argument('--maf', type=float, default=0.05)
    p.add_argument('--max-missing', type=float, default=0.10)
    p.add_argument('--min-count', type=int, default=5)
    p.add_argument('--promoter-length', type=int, default=2000)
    try:
        fetch(p.parse_args())
    except (ValueError, KeyError, requests.RequestException) as e:
        print(f'Collection stopped: {e}. Inspect archived responses; do not substitute demo data.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
