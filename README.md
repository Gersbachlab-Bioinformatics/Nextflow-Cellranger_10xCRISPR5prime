# Nextflow-Cellranger_10xCRISPR5prime
Nextflow Cell Ranger-based pipelines for Perturb-seq experiments.
We have implemented the changes needed to enable using the Cell Ranger package from 10X Genomics to process and create counts tables for transcriptomic gene expression (GEX) and gRNA sequencing data.

## Input files
- gene_exp_guide_lib.csv
  
<img width="353" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/cf02af03-c968-4477-b303-bd668f292f4e">

- guide_feature_reference.csv
  
<img width="593" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/e2fab8e5-5ffc-4a5a-b172-cff3571333e7">

- data
  
<img width="264" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/92551277-7ce3-43d8-8648-eef757d92411">

<img width="263" alt="image" src="https://github.com/Gersbachlab-Bioinformatics/Nextflow-Cellranger_10xCRISPR5prime/assets/104788472/2cd55e59-f3d2-4e74-9630-df8b62402046">

### Pipeline
```
nextflow run /data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/main_cellranger.nf -c /data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/cellranger_perturb.config -with-timeline /data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/cellranger_timeline/cellranger_out -resume -w /data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/cellranger_timeline/cellranger_all
```
