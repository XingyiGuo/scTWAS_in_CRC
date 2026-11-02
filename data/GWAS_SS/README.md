# CRC GWAS summary statistics

The European-ancestry and multi-ancestry colorectal cancer GWAS meta-analyses generated in
this study are deposited in the NHGRI-EBI GWAS Catalog. The summary statistics are not
mirrored in this repository: each meta-analysis is several GB, and the GWAS Catalog is the
citable, versioned source.

| Meta-analysis | Accession |
|---|---|
| CRC, European ancestry | *(GWAS Catalog accession -- to be added on release)* |
| CRC, multi-ancestry | *(GWAS Catalog accession -- to be added on release)* |

Previously published European-ancestry CRC GWAS summary statistics used for comparison are
available under accession **GCST90129505**.

Code that consumes these files -- colocalization, SMR/HEIDI and the association testing in
`code/04_association_testing/` and `code/05_colocalization_smr_heidi/` -- expects the
per-chromosome, rsID-annotated dosage format described in Methods.
