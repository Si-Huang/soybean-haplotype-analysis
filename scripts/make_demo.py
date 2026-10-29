"""Create small synthetic exercises for both analysis branches, never real gene data."""
import argparse
from pathlib import Path
import pandas as pd
from common import empty_output, write_json
from style import POPULATIONS


def make_demo(directory):
    root = empty_output(directory)
    for role in ['target', 'tf']:
        out = root / role
        out.mkdir()
        keys = ['101', '102', '103']
        # First two SNPs occupy one codon: AGG -> TCG, jointly Arg -> Ser.
        alleles = [('A', 'G', 'G'), ('T', 'C', 'G'), ('A', 'G', 'A')]
        haps = pd.DataFrame([{'haplotype_id': f'H{i}', 'source_id': f'demo_{i}',
                              **dict(zip(keys, values))} for i, values in enumerate(alleles, 1)])
        haps.to_csv(out / 'haplotypes.csv', index=False)
        sites = pd.DataFrame({'site_id': keys, 'position': [101, 102, 103],
                              'relative_position': [-300, -299, -298],
                              'cds_position': [4, 5, 6], 'alleles': ['A/T', 'G/C', 'G/A']})
        sites.to_csv(out / 'sites.csv', index=False)
        rows, meta, phenotype = [], [], []
        for p, pop in enumerate(POPULATIONS):
            for h in range(1, 4):
                for j in range([6, 4, 3][(h + p) % 3]):
                    sample = f'DEMO_{role}_{p}_{h}_{j}'
                    rows.extend([{'sample_id': sample, 'haplotype_id': f'H{h}', 'population': pop}] * 2)
                    meta.append({'sample_id': sample, 'longitude': 100 + 5 * p + j, 'latitude': 30 + 2 * h + j / 2})
                    phenotype.extend([{'sample_id': sample, 'trait': 'Synthetic response', 'value': ['R', 'S', 'N'][(j + h) % 3]},
                                      {'sample_id': sample, 'trait': 'Synthetic score', 'value': str((j + p) % 4)}])
        assignments = pd.DataFrame(rows)
        assignments.to_csv(out / 'assignments.csv', index=False)
        assignments.groupby(['haplotype_id', 'population']).size().rename('count').reset_index().to_csv(out / 'counts.csv', index=False)
        pd.DataFrame(meta).to_csv(out / 'metadata.csv', index=False)
        pd.DataFrame(phenotype).to_csv(out / 'phenotypes.csv', index=False)
        if role == 'tf':
            (out / 'reference_cds.fasta').write_text('>SYNTHETIC_TRANSCRIPT\nATGAGGGCCTAA\n', encoding='utf-8')
            pd.DataFrame([{'name': 'Demo region', 'start': 1, 'end': 2, 'source': 'Synthetic illustration only'}]).to_csv(out / 'domains.csv', index=False)
        write_json(out / 'manifest.json', {'gene_id': 'SYNTHETIC_NOT_A_REAL_GENE', 'symbol': 'Synthetic exercise',
                                         'role': role, 'assembly': 'synthetic', 'strand': '+', 'synthetic': True,
                                         'site_columns': keys, 'cds_phase_valid': role == 'tf',
                                         'scope_note': 'Invented practice data; no biological conclusions.'})
    print(f'Synthetic exercises written to {root}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', required=True)
    make_demo(Path(p.parse_args().out))
