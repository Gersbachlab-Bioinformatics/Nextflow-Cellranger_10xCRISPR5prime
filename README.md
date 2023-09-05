# Nextflow-Cellranger_10xCRISPR5prime
Nextflow Cell Ranger-based pipelines for Perturb-seq experiments.

We have implemented the changes needed to enable using the Cell Ranger package from 10X Genomics to process and create counts tables for transcriptomic gene expression (GEX) and gRNA sequencing data.

### Software
Nextflow version 23.04.1

Cellranger-7.1.0

### Input files
- [`gene_exp_guide_lib.csv`](gene_exp_guide_lib.csv)
<img width="500" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/19511894-6fbd-4715-9ea3-dad2cbf7c2a3">


- [`guide_feature_reference.csv`](guide_feature_reference.csv)
<img width="1000" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/e2fab8e5-5ffc-4a5a-b172-cff3571333e7">

- [`data`](data)
<img width="800" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/92551277-7ce3-43d8-8648-eef757d92411">

<img width="800" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/2cd55e59-f3d2-4e74-9630-df8b62402046">
  
### Pipeline
```
nextflow run /data/main_cellranger.nf -c /data/cellranger_perturb.config -with-timeline /data/cellranger_timeline/cellranger_out -resume -w /data/cellranger_timeline/cellranger_all
```
