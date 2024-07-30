## Background

This repositorty contains code for performing bioinformatic analysis of *P. falciparum* next-generation sequencing (NGS) data. The repo contains multiple analysis pipelines. The GeoPrediction and MaRS analyses can be run with amplicon or selective whole genome sequencing (sWGA) data, while the sWGA VCF pipeline should only be run with sWGA data. 

All code in this repository is designed to run on CDC SciComp environment from the CycloDeploy Group Working Area.

Various example commands are provided in each section of this readme to run different types of analyses. These scripts are available in the __qsub script__ and I recommend using the qsub script to run the analysis, as it will point to the correct config file, etc. 

#### I'll need to figure out how to point to some of the reference and output files, since there are already outputs there and some of the folders are really large already

### Remove Human Reads

Human derived reads must be removed from the sWGA data before any analysis is run. Likewise, sWGA data uploaded to NCBI SRA must have human-dreived reads removed. All workflows incorporate human-read removal steps (__make sure__) and human-removed reads will be copied to the xyz directory (__make sure__). You can also run the following script to remove human reads without running any of the full analyses:

`nextflow run main.nf -entry mapHuman_then3D7 -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt`

This will remove human reads and map to the 3D7 reference genome to ascertain depth of coverage accross the *P. falciparum* 3D7 reference genome.
#### Also stuff about genome coverage with R script (update this)

### Geo Prediction

#### Prediction with pfs47 and cpmp

#### Prediction with 221 SNP barcode (TBD)

### Malaria Resistance Surveillance (MaRS)

### Pfal_sWGA_VCF

This repository contains a variant calling workflow (VCF) for selective whole genome amplification (sWGA) data generated for *Plasmodium falciparum*. The pipeline is a nextflow conversion of the the bash pipeline from [Niare et al. 2023](https://link.springer.com/article/10.1186/s12936-023-04632-0) (see their [GitHub page](https://github.com/Karaniare/Optimized_GATK4_pipeline/tree/main)). 

Our code also uses a container created by Niare et al.; however, some adaptations were required to run on CDC's SciComp in nextflow. Therefore, you must run the code in CDC's SciComp environment and this code is not meant to be deployed into other environments. 

Run analysis
    Human read removal will automatically occur in this pipeline

### SNIPPY

#### TBD

### Clustering

#### TBD

### Report Generation

#### TBD

### hrp2/3 deletions

### plasmepsin 2 (copy number)

### pfcoronin

### mdr1 copy number

### moi

### Pfmdt and PfTetQ

#### doxycycline resistance from this [paper](https://wwwnc.cdc.gov/eid/article/21/8/15-0524_article#:~:text=falciparum%20metabolite%20drug%20transporter%20gene,was%20later%20confirmed%20(7))



