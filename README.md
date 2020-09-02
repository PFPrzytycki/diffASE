diffASE User’s Guide
================

This guide provides an overview and example of how to calculate
differential Allele-Specific Expression (ASE) as described in
“[Differential Allele-Specific Expression Uncovers Breast Cancer Genes
Dysregulated by Cis Noncoding
Mutation](https://www.cell.com/cell-systems/fulltext/S2405-4712\(20\)30029-6).”

In order to use the method, simply source “diffASEfunctions.R” in R. The
functions to compute differential ASE depend on the R packages: “metap”,
“limma”, and “biomaRt”

``` r
source("diffASEfunctions.R")
```

The main data required to calculate differential ASE are allele counts
at heterozygous sites in a tumor and matched normal sample. These can be
found using the package “AllelicImbalance” in R with the function
“getAlleleCounts” called on BAM files for each sample. Make sure to
use the same set of heterozygous sites for both the tumor and matched
normal files. Heterozygous sites can be found using
“scanForHeterozygotes” from the same package.

Example allele counts should look as follows, with columns for Entrez
gene id, sample name, genomic location, normal sample primary allele
count, normal sample total read count, tumor sample primary allele
count, tumor sample total read count:

| V1 | V2      | V3             | V4 | V5 | V6 | V7 |
| :- | :------ | :------------- | -: | -: | -: | -: |
| 1  | sample1 | chr19:58346952 | 17 | 25 | 19 | 26 |
| 1  | sample2 | chr19:58346952 | 12 | 19 | 21 | 25 |

Using just this ASE data we can compute baseline differential ASE:

``` r
computeASEbaseline(testASE)
```

| sample  | gene |      pval |       ASE |
| :------ | ---: | --------: | --------: |
| sample1 |    1 | 0.5789254 | 0.0507692 |
| sample2 |    1 | 0.0307450 | 0.2084211 |

Next, to compute the purity adjusted differentual ASE, we need a tumor
purity for each sample. For TCGA data many of these have been
pre-computed in the paper “Systematic pan-cancer analysis of tumour
purity.” These should be stored in a table with the sample name in the
first column and the purity in the second.

| V1      |  V2 |
| :------ | --: |
| sample1 | 0.8 |
| sample2 | 0.9 |

Given that purity data, we can compute diffASE-purity:

``` r
computeASEpurity(testASE, allPurity)
```

| sample  | gene |      pval |       ASE |
| :------ | ---: | --------: | --------: |
| sample1 |    1 | 0.5349559 | 0.0634615 |
| sample2 |    1 | 0.0227733 | 0.2315789 |

Finally, expression adjusted differential ASE can be computed if
expression values for each gene in each sample are known. This requires
two tables: one for cpm counts for the normal sample and one with cpm
counts for the matched tumor sample. These should have the Ensembl name
of the gene as the row name and a column for each sample.

|                 | sample1 | sample2 |
| :-------------- | ------: | ------: |
| ENSG00000121410 |     1.2 |     2.1 |

|                 | sample1 | sample2 |
| :-------------- | ------: | ------: |
| ENSG00000121410 |       2 |     1.9 |

With that tumor purity and expression data, we can comute diffASE-exp:

``` r
computeASEexp(testASE, allPurity, cpmTumor, cpmNormal)
```

| sample  | gene |      pval |       ASE |
| :------ | ---: | --------: | --------: |
| sample1 |    1 | 0.4967081 | 0.0761538 |
| sample2 |    1 | 0.0234966 | 0.2291540 |

For comparison to differential ASE, we can also compute tumor sample ASE
to see that sometimes ASE looks significant without the context of the
matched tumor sample:

``` r
computeTumorSampleASE(testASE)
```

| sample  | gene |      pval |       ASE |
| :------ | ---: | --------: | --------: |
| sample1 |    1 | 0.0186029 | 0.2307692 |
| sample2 |    1 | 0.0006739 | 0.3400000 |
