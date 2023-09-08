# Nextflow-Cellranger_10xCRISPR5prime
Nextflow Cell Ranger-based pipelines for Perturb-seq experiments.

We have implemented the changes needed to enable using the Cell Ranger package from 10X Genomics to process and create counts tables for transcriptomic gene expression (GEX) and gRNA sequencing data.

### Software
Nextflow version 23.04.1

Cellranger-7.1.0

### Schema 
<img width="558" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/ef1f9fdc-ded3-45ef-bf9c-8249af8af863">


### Input files
- [`gene_exp_guide_lib.csv`](gene_exp_guide_lib.csv)
<img width="500" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/19511894-6fbd-4715-9ea3-dad2cbf7c2a3">


- [`guide_feature_reference.csv`](guide_feature_reference.csv)
<img width="1000" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/e2fab8e5-5ffc-4a5a-b172-cff3571333e7">


- [`data`](data)
  https://www.synapse.org/#!Synapse:syn52393737
<img width="800" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/fbbc1f73-3a83-4327-a793-4face5fef821">

<img width="800" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/5f439e71-1338-4ae5-b755-9e58b4d4d826">

  
### Pipeline
```
nextflow run /data/main_cellranger.nf -c /data/cellranger_perturb.config -with-timeline /data/cellranger_timeline/cellranger_out -resume -w /data/cellranger_timeline/cellranger_all
```
