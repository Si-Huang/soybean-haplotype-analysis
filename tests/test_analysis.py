"""Scientific invariants: strand, joint codons, independent units, missingness."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from common import cds_features, promoter_window, reconstruct_proteins
from analyze import analyze, association, materials, phenotype_summary, validate
from make_demo import make_demo
from fetch_soyomics import select_variants


class ScientificChecks(unittest.TestCase):
    def test_unknown_missingness_deferred_to_service_not_dropped(self):
        variants = [{'position': 1, 'maf': .2}, {'position': 2, 'maf': .2, 'miss': .5},
                    {'position': 3, 'maf': .2, 'miss': .01}, {'position': 4, 'maf': .001}]
        result = select_variants(variants, 'tf', {1: 1, 2: 2, 3: 3, 4: 4}, .05, .1)
        self.assertEqual([v['position'] for v in result], [1, 3])
        self.assertNotIn('miss', result[0])

    def test_strand_and_joint_codon(self):
        self.assertEqual(promoter_window(100, 200, '+', 2000), (1, 99))
        self.assertEqual(promoter_window(100, 200, '-', 2000), (201, 2200))
        sites = pd.DataFrame({'site_id': ['a', 'b'], 'cds_position': [4, 5], 'alleles': ['A/T', 'G/C']})
        haps = pd.DataFrame([{'haplotype_id': 'H1', 'a': 'T', 'b': 'C'}])
        proteins, reference = reconstruct_proteins(haps, sites, 'ATGAGGTAA', '+')
        self.assertEqual(reference, 'MR*')
        self.assertEqual(proteins['H1'], 'MS*')
        # Genomic alleles on the negative strand must be complemented.
        negative = sites.assign(alleles=['T/A', 'C/G'])
        proteins, _ = reconstruct_proteins(haps.assign(a='A', b='G'), negative, 'ATGAGGTAA', '-')
        self.assertEqual(proteins['H1'], 'MS*')

    def test_unresolved_allele_and_indel_stop_translation(self):
        for allele, allowed in [('N', 'A/G'), ('AT', 'A/AT')]:
            with self.assertRaises(ValueError):
                reconstruct_proteins(pd.DataFrame([{'haplotype_id': 'H1', 's': allele}]),
                                     pd.DataFrame([{'site_id': 's', 'cds_position': 4, 'alleles': allowed}]), 'ATGAGGTAA', '+')

    def test_cds_phase_across_exons(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'ref.gff'
            path.write_text('c\t.\tCDS\t1\t4\t.\t+\t0\tParent=t\nc\t.\tCDS\t10\t14\t.\t+\t2\tParent=t\n')
            chromosome, strand, positions, valid = cds_features(path, 't')
            self.assertTrue(valid)
            self.assertEqual(positions, [1, 2, 3, 4, 10, 11, 12, 13, 14])

    def test_missing_internal_phase_and_unknown_initial_frame(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'ref.gff'
            path.write_text('c\t.\tCDS\t1\t5\t.\t-\t.\tParent=t\nc\t.\tCDS\t10\t13\t.\t-\t0\tParent=t\n')
            _, _, positions, valid = cds_features(path, 't')
            self.assertTrue(valid)
            self.assertEqual(positions, [13, 12, 11, 10, 5, 4, 3, 2, 1])
            path.write_text('c\t.\tCDS\t1\t9\t.\t+\t.\tParent=t\n')
            with self.assertRaisesRegex(ValueError, 'initial CDS phase'):
                cds_features(path, 't')
            path.write_text('c\t.\tCDS\t1\t4\t.\t+\t0\tParent=t\nc\t.\tCDS\t10\t14\t.\t+\t1\tParent=t\n')
            self.assertFalse(cds_features(path, 't')[3])

    def test_materials_and_fieldwise_phenotype_missingness(self):
        assignments = pd.DataFrame([['a', 'H1', 'Soja'], ['a', 'H1', 'Soja'],
                                    ['b', 'H1', 'Soja'], ['b', 'H2', 'Soja'], ['c', 'H2', 'Landrace']],
                                   columns=['sample_id', 'haplotype_id', 'population'])
        unique, excluded = materials(assignments)
        self.assertEqual(set(unique.sample_id), {'a', 'c'})
        self.assertEqual(excluded.sample_id.tolist(), ['b'])
        records = pd.DataFrame([['a', 'response', 'N'], ['a', 'response', 'N'], ['a', 'score', '3.0'],
                                ['a', 'score', '3'], ['c', 'response', 'R'], ['c', 'response', 'S'],
                                ['c', 'score', 'NA']], columns=['sample_id', 'trait', 'value'])
        summary, conflicts, matched = phenotype_summary(unique, records, ['H1', 'H2'])
        self.assertEqual(len(matched), 2)
        self.assertEqual(conflicts.to_dict('records'), [{'sample_id': 'c', 'trait': 'response'}])
        self.assertEqual(matched[matched.trait == 'score'].value.tolist(), ['3'])
        self.assertTrue(summary[summary.haplotype_id == 'H2'].percent.isna().all())
        self.assertEqual(summary[(summary.haplotype_id == 'H1') & (summary.category == 'N')]['count'].iloc[0], 1)

    def test_sparse_reproducible_association(self):
        unique = pd.DataFrame({'sample_id': list('abcdef'), 'haplotype_id': ['H1', 'H1', 'H2', 'H2', 'H3', 'H3'],
                               'population': ['Soja', 'Landrace', 'Soja', 'Improved', 'Landrace', 'Improved']})
        a, b = association(unique, 'CDS'), association(unique, 'CDS')
        self.assertEqual(a['p'], b['p'])
        self.assertEqual(a['n'], 6)
        self.assertIn('Monte Carlo', a['method'])

    def test_copy_count_mismatch_is_rejected(self):
        haps = pd.DataFrame([{'haplotype_id': 'H1', 's': 'A'}])
        sites = pd.DataFrame([{'site_id': 's', 'position': 1, 'alleles': 'A/G'}])
        counts = pd.DataFrame({'haplotype_id': ['H1'] * 3, 'population': ['Soja', 'Landrace', 'Improved'], 'count': [2, 0, 0]})
        assignments = pd.DataFrame([{'sample_id': 'a', 'haplotype_id': 'H1', 'population': 'Soja'}])
        with self.assertRaisesRegex(ValueError, 'do not reproduce'):
            validate(haps, sites, counts, assignments)

    def test_no_variants_does_not_invent_cohort(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / 'input'
            source.mkdir()
            (source / 'manifest.json').write_text(json.dumps({'role': 'target', 'strand': '+', 'assembly': 'synthetic'}))
            pd.DataFrame(columns=['haplotype_id']).to_csv(source / 'haplotypes.csv', index=False)
            pd.DataFrame(columns=['site_id', 'position', 'alleles']).to_csv(source / 'sites.csv', index=False)
            pd.DataFrame(columns=['haplotype_id', 'population', 'count']).to_csv(source / 'counts.csv', index=False)
            analyze(argparse.Namespace(input=str(source), out=str(root / 'out'), phenotypes=None, domains=None))
            audit = json.loads((root / 'out' / 'run_audit.json').read_text(encoding='utf-8'))
            self.assertIn('no retained', audit['modules']['core'])
            self.assertFalse((root / 'out' / 'haplotype_matrix.csv').exists())

    def test_both_branches_end_to_end(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            make_demo(root / 'inputs')
            for role in ['target', 'tf']:
                source, out = root / 'inputs' / role, root / role
                analyze(argparse.Namespace(input=str(source), out=str(out), phenotypes=str(source / 'phenotypes.csv'),
                                           domains=str(source / 'domains.csv') if role == 'tf' else None))
                audit = json.loads((out / 'run_audit.json').read_text(encoding='utf-8'))
                self.assertEqual(audit['modules']['core'], 'ok')
                self.assertTrue((out / 'figures' / 'haplotype_matrix.svg').exists())
                if role == 'tf':
                    mapping = pd.read_csv(out / 'cds_to_protein.csv')
                    self.assertEqual(mapping.protein_type.nunique(), 2)  # synonymous H1/H3 collapse
                    self.assertTrue(audit['modules']['protein'].startswith('ok:'))
                frequencies = pd.read_csv(out / 'haplotype_frequencies.csv')
                self.assertTrue((frequencies.groupby('population').pct_within_population.sum().round(6) == 100).all())
                with self.assertRaises(ValueError):
                    analyze(argparse.Namespace(input=str(source), out=str(out), phenotypes=None, domains=None))


if __name__ == '__main__':
    unittest.main()
