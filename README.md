# ***Dopamine release from Parkinson’s patient-derived neurons is disrupted due to impaired synaptic vesicle loading***

**This repository contains the scripts and metadata necessary to replicate the transcriptomics analyses included in the manuscript.**


---


## Data / code availability

- Raw (FASTQ) and processed (Salmon quant.sf) bulk RNA-seq data: available from the Gene Expression Omnibus (GEO) under accession number GSE353879.
- Metadata: included in this GitHub repository ('benvallin/Cramb_2026'), within folder 'metadata'.
- Scripts: included in this GitHub repository ('benvallin/Cramb_2026'), within folder 'scripts'.
- Software details: included in this GitHub repository ('benvallin/Cramb_2026'), within folder 'software_details'.
- Pseudobulk-level differential gene expression analysis of publicly available snRNA-seq datasets: scripts and results available upon request.


---


## Code execution steps

**RWM dataset**
1. rwm/make_feature_metadata.R
2. rwm/dgea.R
3. rwm/slc18a2_expression.R

**Meta-analysis of publicly available snRNA-seq datasets**
1. public/random_effects_metaanalysis.R


---


## Requirements

- R package `rnaseqtools` available from [GitHub](https://github.com/benvallin/rnaseqtools) with:

```r
pak::pak("benvallin/rnaseqtools")
```

- R package `paramstream` available from [GitHub](https://github.com/benvallin/paramstream) with:

```r
pak::pak("benvallin/paramstream")
```


---