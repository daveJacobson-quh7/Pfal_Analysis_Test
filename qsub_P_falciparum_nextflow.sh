#!/bin/bash -l
#$ -e pfalciparum_nextflow.err
#$ -o pfalciparum_nextflow.out
#$ -N pfalciparum_nextflow
#$ -q cyclospora.q
#$ -cwd
#$ -pe smp 8


myDir=/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/P_falciparum_fullAnalysis
export LC_ALL=C
#export JAVA_CMD=/usr/bin/java
# export JAVA_CMD=/apps/x86_64/java/openjdk1.8.0_41/bin/java

module load nextflow

myDate=$(date +"%Y-%m-%d_%H-%M-%S")

#Modify below this line

#Perform read QC only
# nextflow run pfal_geoMain.nf -entry qcReads -profile conda -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt -resume

#Run the complete pipeline - using both cpmp and pfs47

# nextflow run main.nf -entry P_falciparum_GeoClassifier_complete_pipeline -profile singularity  -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt 

# nextflow run main.nf -entry P_falciparum_GeoClassifier_useCleanReads -profile singularity  -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt 

# nextflow run main.nf -entry P_falciparum_MaRS -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt 

# nextflow run main.nf -entry P_falciparum_MaRS_GeoPrediction -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt --resume

# nextflow run main.nf -entry P_falciparum_GeoClassifier_useExistingVCF -profile singularity  -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt --metadata $myDir/sWGS_metadata.xlsx

# nextflow run main.nf -entry mapHuman_then3D7 -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt

# nextflow run main.nf -entry just3D7 -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt

# nextflow run main.nf -entry myDepth -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt -resume

# nextflow run main.nf -entry optimizedGATK4_VCF -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt -resume


