// Declare syntax version
nextflow.enable.dsl=2
// Script parameters
//params.GTF_GZ_LINK = 'http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz'
//params.CELLRANGER_REF = '/data/gersbachlab/Ruhi/pipeline_perturbseq_like/refdata-gex-GRCh38-2020-A'
//params.TRANSCRIPTOME_REFERENCE = "human"
//params.KALLISTO_BIN = '/data/reddylab/Ruhi/miniconda3/envs/perturbseq_like_pipeline/bin/kallisto'
//params.GENOME = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
//params.THREADS = 15
//params.DISTANCE_NEIGHBORS = 1000000
//params.IN_TRANS = "FALSE"
//params.FASTQ_FILES_TRANSCRIPTS = ['/data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L001_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L001_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L002_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L002_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L003_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L003_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L004_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRa_D1_S4_L004_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L001_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L001_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L002_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L002_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L003_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L003_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L004_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/GEX_CRISPRi_D1_S1_L004_R2_001.fastq.gz']
//params.FASTQ_NAMES_TRANSCRIPTS = ['S1_L1']
//params.FASTQ_FILES_GUIDES = ['/data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L001_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L001_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L002_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L002_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L003_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L003_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L004_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRa_D1_S10_L004_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L001_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L001_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L002_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L002_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L003_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L003_R2_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L004_R1_001.fastq.gz /data/reddylab/Alex/collab/20220523_Sean/data/sc_rna_seq/fastq/test/gRNA_CRISPRi_D1_S7_L004_R2_001.fastq.gz']
//params.FASTQ_NAMES_GUIDES = ['S7_L1']
//params.WHITELIST= '/data/gersbachlab/Ruhi/Nextflow/737K-august-2016.txt'


workflow {
             
    cellranger_out = cellranger (
                 Channel.from(params.FASTQ_NAMES_TRANSCRIPTS),
                 Channel.from(params.FASTQ_FILES_TRANSCRIPTS),
                 Channel.from(params.FASTQ_NAMES_GUIDES),
                 Channel.from(params.FASTQ_FILES_GUIDES),
                 Channel.from(params.CELLRANGER_REF),
                 Channel.of(params.THREADS).collect(),
                 Channel.of(params.WHITELIST).collect(),
                )
    
}


process cellranger {
        input:
        tuple val(rna_sample_name)
        tuple val(rna_string_fastqz)
        tuple val(guide_sample_name)
        tuple val(guide_string_fastqz)
        tuple val(cellranger_ref)
        tuple val(threads)
        tuple val(whitelist)
        output:
        path ('outs/filtered_feature_bc_matrix'),  emit: mtx_file

        script:
    
        """
        
        chmod 700 /data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/generating_cell_range_inputs.py
        /data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/generating_cell_range_inputs.py --rna_sample_name $rna_sample_name --rna_string_fastqz "$rna_string_fastqz" --guide_sample_name $guide_sample_name --guide_string_fastqz "$guide_string_fastqz" --cellranger_ref $cellranger_ref --threads 15 --whitelist $whitelist
        #The script generating_cell_range_inputs.py should create the library.csv and feature_ref.csv 
        #remove the comment to run cellranger
        /data/gersbachlab/Ruhi/pipeline_perturbseq_like/cellranger-7.1.0/bin/cellranger count --id=20210523_Sean --libraries=/data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/gene_exp_guide_lib.test.csv --transcriptome=$cellranger_ref --feature-ref=/data/gersbachlab/Ruhi/Nextflow/IGVF_Nextflow/guide_feature_reference.csv
        
        #Artificially creating the output. Remove it when cellranger starts running.
        
        mkdir outs
        mkdir outs/filtered_feature_bc_matrix/
        mkdir 'outs/filtered_feature_bc_matrix/matrix.mtx.gz'
        
      
        
        """
} 



