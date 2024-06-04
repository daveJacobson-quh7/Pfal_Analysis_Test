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

# nextflow run mainCombo.nf -entry P_falciparum_GeoClassifier_complete_pipeline -profile singularity  -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt 

# nextflow run mainCombo.nf -entry P_falciparum_GeoClassifier_useCleanReads -profile singularity  -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt 

# nextflow run mainCombo.nf -entry P_falciparum_MaRS -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt 

nextflow run mainCombo.nf -entry P_falciparum_MaRS_GeoPrediction -profile singularity -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt --resume

# nextflow run pfal_geoMain.nf  -entry P_falciparum_classifier_complete_pipeline -profile conda -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt --balkClassifier /scicomp/home-pure/yyr4/Desktop/balk_classifier/lkc-bisnp/lkc_bisnp-run --balkRegion_JLIB  /scicomp/home-pure/yyr4/Desktop/balk_classifier/lkc-bisnp/pf_combinedDhruvi_trainRegion.joblib.gz --balkCountry_JLIB  /scicomp/home-pure/yyr4/Desktop/balk_classifier/lkc-bisnp/pf_combinedDhruvi_trainCoumtry.joblib.gz --metadata $myDir/DMS18_metadata_Geo.xlsx -resume 
#Run the complete pipeline - use only pfs47
# nextflow run pfal_geoMain.nf  -entry P_falciparum_classifier_complete_pipeline -profile singularity -resume -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt --balkContinent_JLIB $myDir/REFERENCES/pfs47_1_continent_classifier.joblib.gz --balkCountry_JLIB $myDir/REFERENCES/pfs47_1_country_classifier.joblib.gz --balkRegion_JLIB $myDir/REFERENCES/pfs47_1_region_classifier.joblib.gz   --ampliseqRef $myDir/REFERENCES/Pfs47_ref.fasta --panelRefBarcode $myDir/REFERENCES/pfs47_barcode_allSamples.csv --snpPositions $myDir/REFERENCES/pfs47_snpPositions.csv --balkKeyFile $myDir/REFERENCES/balkKey_pfs47.txt --ampliconList $myDir/REFERENCES/pfs47_amplicon.txt

# nextflow run mainCombo.nf -entry P_falciparum_GeoClassifier_useExistingVCF -profile singularity  -with-report nextflow_logFiles/$myDate\_report.html -with-trace nextflow_logFiles/$myDate\_trace.txt --metadata $myDir/sWGS_metadata.xlsx


