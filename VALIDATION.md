# Validation

Checked on 2026-10-29 with Python 3.14 on Windows.

- Skill frontmatter validated with the Codex skill-creator validator.
- Unit/integration checks cover strand orientation, CDS exon phase, joint codon changes, synonymous protein collapse, missing phenotype values, field-specific conflicts, independent-material counting, summary/assignment mismatch, empty variant sets, sparse-table statistics, and both complete synthetic workflows.
- Live SoyOmics retrieval and analysis exercised both branches: LOX1w (`SoyZH13_13G319900`, promoter) and SVPj (`SoyZH13_15G055900.m1`, CDS/protein). SVPj reference phase/length and SNP protein reconstruction passed. LOX1w includes multi-base alleles, used to check matrix column sizing.
- Matrix, population, protein-position, network, and synthetic phenotype figures inspected locally. Real material coordinates and phenotype records were not available in these live bundles; those modules correctly reported unavailable inputs. Synthetic fixtures exercise these optional modules.
- GitHub Actions runs the offline scientific checks on Python 3.12. Live service availability is separate from the offline suite.

No real material-level data or local analysis outputs are included in this repository. Example images are generated entirely from synthetic exercises.
