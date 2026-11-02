# CRC GWAS summary statistics

The European-ancestry and multi-ancestry colorectal cancer GWAS meta-analyses generated in
this study are deposited in the NHGRI-EBI GWAS Catalog, and are also included in the Zenodo
record for this release (see the DOI in the top-level `README.md`, under
`CRC_GWAS_SS/`). They are not mirrored in this repository: together they are over 800 MB,
and git is the wrong place for them.

The per-chromosome files are `final_gc_maf0.01_CRC_{EUR,Multi}_META_1.chr<N>.chrpos.rsid.assoc.dosage.gz`
(22 autosomes per ancestry), rsID-annotated, with the lead-variant table
`CRC_238signals_plus_additional_EUR_MULTI_loci.csv` alongside.

| Meta-analysis | Accession |
|---|---|
| CRC, European ancestry | *(GWAS Catalog accession -- to be added on release)* |
| CRC, multi-ancestry | *(GWAS Catalog accession -- to be added on release)* |

Previously published European-ancestry CRC GWAS summary statistics used for comparison are
available under accession **GCST90129505**.

Code that consumes these files -- colocalization, SMR/HEIDI and the association testing in
`code/04_association_testing/` and `code/05_colocalization_smr_heidi/` -- expects the
per-chromosome, rsID-annotated dosage format described in Methods.
