

process trainBALK {
  label 'balkClassifier'

  input:
  
  output:
  path "*Region*.joblib.gz"
  path "*Country*.joblib.gz"

  script:
  """
  ${params.balkClassifier} train --posfile ${params.snpPositions} --method balk -o geoCombined_finalCountry_classifier.joblib.gz $baseDir/REFERENCES/geoCombined_trainCountryRenameFinal_1.csv
  ${params.balkClassifier} train --posfile ${params.snpPositions} --method balk -o geoCombined_finalRegion_classifier.joblib.gz $baseDir/REFERENCES/geoCombined_trainRegionRenameFinal_1.csv


  """
}




process fastqc_beforeTrim {

  // label 'fastqc_container'

  // label 'fastqcEnv'

  label 'bioinformaticProcessing'

  publishDir "${params.qcFolder}", pattern: "*qualityCheck.txt", mode: "copy"



  input:
  tuple val(sample_id), path(sample_files)

  output:
  path "*qualityCheck.txt"

  shell:

  '''

  mkdir fastqc_out

  fastqc !{sample_files} -o fastqc_out

  for x in fastqc_out/*.zip; do unzip "$x"; done

  echo !{sample_files[0].simpleName}_fastqc > myName



  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files[0].simpleName}_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_id}.clean1_qualityCheck.txt   
  grep "Total Sequences" !{sample_files[0].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean1_qualityCheck.txt
  grep "Sequence length" !{sample_files[0].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean1_qualityCheck.txt

  

  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files[1].simpleName}_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_id}.clean2_qualityCheck.txt  
  grep "Total Sequences" !{sample_files[1].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean2_qualityCheck.txt
  grep "Sequence length" !{sample_files[1].simpleName}_fastqc/fastqc_data.txt >> !{sample_id}.clean2_qualityCheck.txt



  '''



}



process mapPhix {

  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

 // publishDir "${qcFolder}", pattern: "*phix.sam", mode: "copy"
  publishDir "${params.qcFolder}", pattern: "*phix.out", mode: "copy"


  input:
  tuple val(sample_id), path(sample_files)

  output:
  path "*phix.out"
  //path "*phix.sam"



  script:

  """

  bowtie2 -p 5 -x $projectDir/REFERENCES/illumina_phiX.bt_index -1 ${sample_files[0]} -2 ${sample_files[1]} -S ${sample_id}.phix.sam &> ${sample_id}.phix.out

  """
}



process trimReads {
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

  input:
  tuple val(sample_id), path(sample_files)


  output:
  tuple val(sample_id), path("*.bbduk.fastq.gz"), emit: trimmedFASTQ
  tuple val(sample_id), path("${sample_id}.stats.txt"), emit: Trimmed_stats


  script:


   """

   bbduk.sh  -Xmx${params.bbmergeRAM} ktrim=r k=27 mink=11 hdist=1 edist=0 ref=${params.adapters} qtrim=rl trimq=20 minlength=50 trimbyoverlap=t minoverlap=24 qin=33 in=${sample_files[0]} in2=${sample_files[1]} out=${sample_id}.R1.bbduk.fastq out2=${sample_id}.R2.bbduk.fastq stats=${sample_id}.stats.txt

   gzip *.fastq
  
  """
}


// bbduk.sh  -Xmx1g  ktrimright=t k=27 hdist=1 edist=0 ref=${adapters} \\
//      qtrim=rl trimq=35 minlength=100 trimbyoverlap=t minoverlap=24 qin=33 in=${reads[0]} in2=${reads[1]} \\
//      out=${sample_id}.R1.bbduk.fastq out2=${sample_id}.R2.bbduk.fastq stats=${sample_id}.stats.txt

process fastqc_afterTrim {

  // label 'fastqc_container'
  label 'bioinformaticProcessing'

  // label 'fastqcEnv'

  publishDir "${params.qcFolder}", pattern: "*qualityCheck.txt", mode: "copy"



  input:
  tuple val(sample_id), path(sample_files)

  output:
  path "*qualityCheck.txt"

  shell:

  '''

  mkdir fastqc_out

  fastqc !{sample_files} -o fastqc_out

  for x in fastqc_out/*.zip; do unzip "$x"; done

  echo !{sample_files[0].simpleName}_fastqc > myName



  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files[0].simpleName}.R1.bbduk_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_id}.clean1_postTrim_qualityCheck.txt   
  grep "Total Sequences" !{sample_files[0].simpleName}.R1.bbduk_fastqc/fastqc_data.txt >> !{sample_id}.clean1_postTrim_qualityCheck.txt
  grep "Sequence length" !{sample_files[0].simpleName}.R1.bbduk_fastqc/fastqc_data.txt >> !{sample_id}.clean1_postTrim_qualityCheck.txt

  

  echo "Average Quality" `sed -n '/^#Base.*Mean/,/^>>END_MODULE/p' !{sample_files[1].simpleName}_fastqc/fastqc_data.txt | grep -v "^>" | grep -v "^#B" |  awk '{print$2}' | awk '{s+=$1}END{print s/NR}'` > !{sample_id}.clean2_postTrim_qualityCheck.txt  
  grep "Total Sequences" !{sample_files[1].simpleName}.R2.bbduk_fastqc/fastqc_data.txt >> !{sample_id}.clean2_postTrim_qualityCheck.txt
  grep "Sequence length" !{sample_files[1].simpleName}.R2.bbduk_fastqc/fastqc_data.txt >> !{sample_id}.clean2_postTrim_qualityCheck.txt



  '''



}



process humanMap{

  // label 'mappingSoftware'
  label 'bioinformaticProcessing'


  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*.bothReadsUnmapped_sorted.bam", emit: sortedNoHumanBAM


  script:
  """
  bowtie2 -p 8 -x ${params.humanRef} -1 ${sample_files[0]} -2 ${sample_files[1]}   -S ${sampleID}_mapped_and_unmapped.sam

  echo 'after human mapping'

  samtools view -bS ${sampleID}_mapped_and_unmapped.sam > ${sampleID}_mapped_and_unmapped.bam

  samtools view -b -f 12 -F 256 ${sampleID}_mapped_and_unmapped.bam  > ${sampleID}_bothReadsUnmapped.bam 

  samtools sort -n -m 5G -@ 2 ${sampleID}_bothReadsUnmapped.bam -o ${sampleID}.bothReadsUnmapped_sorted.bam
  """
}


process noHumanReads{
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

  publishDir "${params.cleanReads}", pattern: "*.host_removed*.fastq.gz", mode: "copy"
  
  input:
  path noHumanBAM

  output:
  tuple val(noHumanBAM.simpleName), path("*.host_removed*.fastq.gz"), emit: humanRemoveFASTQ
  path("*.host_removed*.fastq.gz")

  script:
  """
 
  samtools fastq -@ 8 ${noHumanBAM} -1 ${noHumanBAM.simpleName}.host_removed_R1.fastq.gz -2 ${noHumanBAM.simpleName}.host_removed_R2.fastq.gz -0 ${noHumanBAM.simpleName}_others.fastq.gz -s ${noHumanBAM.simpleName}_singletons.fastq.gz -n
  
  """
}


process indexRefs{
  label 'bioinformaticProcessing'

  script:
  """

  samtools faidx ${params.ampliseqRef}
  gatk CreateSequenceDictionary -R ${params.ampliseqRef}
  """
}



process mapToRefs {

  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*ampliseqMap.sam", emit: mappedSAM

  script:
  """
  bwa mem -t 16 -M -R  "@RG\\tID:'${sampleID}'\\tLB:'${sampleID}'\\tPL:illumina\\tSM:'${sampleID}'\\tPU:'${sampleID}'" ${params.ampliseqRef} ${sampleID}.host_removed_R1.fastq.gz ${sampleID}.host_removed_R2.fastq.gz > ${sampleID}.ampliseqMap.sam

  """
}

process mapToRefs_bowtie2 {
  // errorStrategy 'ignore'
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*mapped_sorted_merged.bam", emit: bowtieMapped

  shell:
  '''
  cat !{params.ampliconList} | while read line; do bowtie2 -x !{params.REFERENCES}/$line\\_bt2 -1 !{sample_files[0].simpleName}.host_removed_R1.fastq.gz -2 !{sample_files[0].simpleName}.host_removed_R2.fastq.gz --local  --rg-id !{sample_files[0].simpleName} --rg SM:!{sample_files[0].simpleName} --rg LB:!{sample_files[0].simpleName} --rg PU:!{sample_files[0].simpleName} --rg PL:ILLUMINA !{sample_files[0].simpleName}_$line.sam ;done
  
  cat !{params.ampliconList} | while read line; do samtools view -b -F 4 !{sample_files[0].simpleName}_$line.sam > !{sample_files[0].simpleName}_$line\\_mapped.bam ; done

  
  samtools merge !{sample_files[0].simpleName}_mapped_merged.bam !{sample_files[0].simpleName}_*mapped.bam 

  samtools sort !{sample_files[0].simpleName}_mapped_merged.bam > !{sample_files[0].simpleName}.mapped_sorted_merged.bam

  '''
}



process cleanSAM_bowtie2{
  label 'bioinformaticProcessing'

  input:
  path mapBowtieSAM

  output:
  path "*mapped_sorted_merged_clean.bam", emit: cleanedSAM

  script:
  """
  gatk  CleanSam -R ${params.ampliseqRef} -I ${mapBowtieSAM.simpleName}.mapped_sorted_merged.bam -O ${mapBowtieSAM.simpleName}.mapped_sorted_merged_clean.bam
  """

}

process indexBAM_bowtie2{
  // label 'variantCalling'
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'
  input:
  path markedBAM
  output:
  tuple val(markedBAM.simpleName), path("*.mapped_sorted_merged_clean_index*"), emit: indexOut
  script:
  """
  cp ${markedBAM} ${markedBAM.simpleName}.mapped_sorted_merged_clean_index.bam
  samtools index ${markedBAM.simpleName}.mapped_sorted_merged_clean_index.bam
  """
  
}

process gatkHapCaller_bowtie2{
  //label 'variantCalling'
  label 'bioinformaticProcessing'
  publishDir "${params.variantFolder}", pattern: "*gatk_bowtie2.g.vcf", mode: "copy"
  input:
  tuple val(sampleID), path(markedDupBAM)
  output:
  path "*gatk_bowtie2.g.vcf", emit: gatkOut
  script:
  """
  gatk  HaplotypeCaller -R ${params.ampliseqRef} -I ${sampleID}.mapped_sorted_merged_clean_index.bam --native-pair-hmm-threads 8  -ERC GVCF  -O ${sampleID}.ampliseq_gatk_bowtie2.g.vcf 
  """

}

process cleanSAM {

 // label 'variantCalling'
 label 'bioinformaticProcessing'

  input:
  path mapSAM

  output:

  path "*ampliseq.clean.bam", emit: cleanedSAM

  script:
  """


  gatk  SamFormatConverter -R  ${params.ampliseqRef} -I ${mapSAM} -O ${mapSAM.simpleName}.ampliseqMap.bam
  gatk  CleanSam -R ${params.ampliseqRef} -I ${mapSAM.simpleName}.ampliseqMap.bam -O ${mapSAM.simpleName}.ampliseq.clean.bam
  """
}

process markDupes {
  //label 'variantCalling'
  label 'bioinformaticProcessing'

  input:
  path cleanBAM

  output:
  path "*ampliseq.sorted.dup.bam", emit: singleFile

  // tuple val(cleanBAM.simpleName), path("*.sorted.dup*"), emit: secondOut

  script:
  """
  gatk MarkDuplicatesSpark -I ${cleanBAM} -O ${cleanBAM.simpleName}.ampliseq.sorted.dup.bam 

  """
}


process indexBAM {
  // label 'variantCalling'
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

  input:
  path markedBAM

  output:
  tuple val(markedBAM.simpleName), path("*.sorted.dup.index*"), emit: indexOut

  script:
  """
  cp ${markedBAM} ${markedBAM.simpleName}.ampliseq.sorted.dup.index.bam
  samtools index ${markedBAM.simpleName}.ampliseq.sorted.dup.index.bam
  """
}

process gatkHapCaller {

  //label 'variantCalling'
  label 'bioinformaticProcessing'
  publishDir "${params.variantFolder}", pattern: "*gatk.g.vcf", mode: "copy"

  input:
  tuple val(sampleID), path(markedDupBAM)


  output:
  path "*ampliseq_gatk.g.vcf", emit: gatkOut

  script:
  """
  gatk  HaplotypeCaller -R ${params.ampliseqRef} -I ${sampleID}.ampliseq.sorted.dup.index.bam --native-pair-hmm-threads 8  -ERC GVCF  -O ${sampleID}.ampliseq_gatk.g.vcf 
  """

}



process gatherGVCFs{
 // label 'variantCalling'
 label 'bioinformaticProcessing'

  input:
  path vcfList

  output:
  tuple val("gatk"), path("gatk.g.vc*"), emit: combinedGVCFs

  script:
  """
  echo "${vcfList.join('\n')}" > vcf.list

  gatk CombineGVCFs -R ${params.ampliseqRef} --variant vcf.list -O gatk.g.vcf
  """
}

process genotypeGATK {
 // label 'variantCalling'
 label 'bioinformaticProcessing'

  input:
  tuple val("gatk"), path(gvcf)

  output:
  path "gatk_allLoci.vcf", emit: gatkVCF

  script:
  """
  gatk GenotypeGVCFs -R ${params.ampliseqRef}  -V gatk.g.vcf --include-non-variant-sites -O gatk_allLoci.vcf

  """

}

process variantsToTable_gatk {
 // label 'variantCalling'
  // publishDir "${params.variantFolder}", pattern: "*gatk_variantTable.table", mode: "copy"
  label 'bioinformaticProcessing'
  input:
  path gatkVCF

  output:
  path "gatk_variantTable.table", emit: gatkTable

  script:
  """
  gatk VariantsToTable -V ${gatkVCF} -F CHROM -F POS -GF GT -GF DP --show-filtered -O gatk_variantTable.table
  """

}

//Fix this to have the ref files better set up
process gatkBarcode {
  publishDir "${params.variantFolder}", pattern: "*gatk_pfal_ampliseq_barcode.csv", mode: "copy"
  publishDir "${params.variantFolder}", pattern: "*.table", mode: "copy"

  // label 'pythonEnv'
  label 'bioinformaticProcessing'

  input:
  path gatkTable
  path vcfList

  output:
  path "*gatk_pfal_ampliseq_barcode.csv", emit:barcodeOut
  path "*.table"

  script:
  """
  echo "${vcfList.join('\n')}" > vcf.list
  sed 's/.ampliseq_gatk_bowtie2.g.vcf//g' vcf.list > sampleList.txt

  python3 $baseDir/python_R_scripts/snpTable_toBarcode_Pf.py --balkKey ${params.balkKeyFile} --snpFile ${gatkTable} --refBarcode ${params.panelRefBarcode} --sampleList sampleList.txt --outfile gatk_pfal_ampliseq_barcode.csv 

  """
}

//Could also include steps for creating the classifier as a separte workflow
//Probably also add a timestamp or something uniq to prediction name - maybe same thing with barcodes above
process predictGeo_gatk {
  label 'balkClassifier'

  input:
  path sampleBarcode

  output:
  path "samples_gatk_region_predicted.out", emit: gatkRegion
  path "samples_gatk_country_predicted.out", emit: gatkCountry
  path "samples_gatk_continent_predicted.out", emit: gatkContinent

  script:
  """

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_gatk_region_predicted.out --jlibfile  ${params.balkRegion_JLIB}

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_gatk_country_predicted.out --jlibfile  ${params.balkCountry_JLIB}

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_gatk_continent_predicted.out --jlibfile  ${params.balkContinent_JLIB}
  """
}



process bcftoolsMpileup {
  // label 'variantCalling'
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'

  input:
  tuple val(sampleID), path(markedDupBAM)

  output:
  path "*mpileup", emit: pileupOut

  script:
  """
  bcftools mpileup -O b -o ${sampleID}.mpileup --annotate FORMAT/DP  -f ${params.ampliseqRef} ${sampleID}.mapped_sorted_merged_clean_index.bam
  """
}

process bcftoolsCall {
  // label 'mappingSoftware'
  // label 'variantCalling'
  label 'bioinformaticProcessing'

  publishDir "${params.variantFolder}", pattern: "*.bcf.vcf", mode: "copy"

  input:
  path mpileupFile

  output:
  path "*bcf.vcf", emit: callOut

  script:
  """
  bcftools call -m  -o ${mpileupFile.simpleName}.bcf.vcf ${mpileupFile}
  """
}


process bcftoolsFilter_Index {
  // label 'variantCalling'
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'
  input:
  path bcfVCF

  output:
  // tuple val(bcfVCF.simpleName), path("*ampliseq_bcftools.vcf.gz"), emit: sampleList
  // tuple val(bcfVCF.simpleName), path("*ampliseq_bcftools.vcf.*"), emit: vcfAndIndex
  path "*ampliseq_bcftools.vcf.gz.csi", emit: vcfIndex
  path "*ampliseq_bcftools.vcf.gz", emit: vcfOnly

  script:
  """
  bcftools filter -O z ${bcfVCF} -o ${bcfVCF.simpleName}.ampliseq_bcftools.vcf.gz
  bcftools index ${bcfVCF.simpleName}.ampliseq_bcftools.vcf.gz
  """

}


process mergeBCFtools{
  // label 'variantCalling'
  // label 'mappingSoftware'
  label 'bioinformaticProcessing'
  input:
  // tuple val(sampleID), path(vcfFiles)
  path vcfsIn
  path indexIn

  output:
  path "bcftools_merge.vcf", emit: bcftoolsMerged

  script:
  """

  echo "${vcfsIn.join('\n')}" > vcf.list


  bcftools merge -l vcf.list -o bcftools_merge.vcf
  """
}

process variantsToTable_bcftools {
  //label 'variantCalling'
  label 'bioinformaticProcessing'

  input:
  path bcfVCF

  output:
  path "bcftools_variantTable.table", emit: bcftoolsTable

  script:
  """
  gatk VariantsToTable -V ${bcfVCF}  -F CHROM -F POS -GF GT -GF DP -O bcftools_variantTable.table
  """

}

//Fix this to have the ref files better set up
process bcftoolsBarcode {
  publishDir "${params.variantFolder}", pattern: "*bcftools_pvivax_ampliseq_barcode.csv", mode: "copy"
  publishDir "${params.variantFolder}", pattern: "*.table", mode: "copy"
  // label 'pythonEnv'
  label 'bioinformaticProcessing'

  input:
  path bcfTable
  path vcfsIn

  output:
  path "*bcftools_pvivax_ampliseq_barcode.csv", emit:barcodeOut
  path "*.table"
  
  script:
  """
  echo "${vcfsIn.join('\n')}" > vcf.list
  sed 's/.ampliseq_bcftools.vcf.gz//g' vcf.list > sampleList.txt

  python3 $baseDir/python_R_scripts/snpTable_toBarcode_Pf.py --balkKey ${params.balkKeyFile} --snpFile ${bcfTable} --refBarcode ${params.panelRefBarcode} --sampleList sampleList.txt --outfile bcftools_pvivax_ampliseq_barcode.csv 

  """
}


process predictGeo_bcftools {
  label 'balkClassifier'

  input:
  path sampleBarcode

  output:
  path "samples_bcftools_region_predicted.out", emit: bcfRegion
  path "samples_bcftools_country_predicted.out", emit: bcfCountry
  path "samples_bcftools_continent_predicted.out", emit: bcfContinent

  script:
  """

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_bcftools_region_predicted.out --jlibfile  ${params.balkRegion_JLIB}

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_bcftools_country_predicted.out --jlibfile  ${params.balkCountry_JLIB}

  ${params.balkClassifier} predict ${sampleBarcode} -o samples_bcftools_continent_predicted.out --jlibfile  ${params.balkContinent_JLIB}
  """
}

//Order is important here
//gatk country, gatk region, bcftools country, bcftools region
//I don't know why, but i had to give the conda path for Rscript, otherwise it would use base R and not R with packages i installed
process mergePredictions {
  // label 'rMerging'
  label 'bioinformaticProcessing'
  publishDir "${params.predictionFolder}", pattern: "*prediction_withMetadata.txt", mode: "copy"

  input:
  path gatkCountry
  path gatkRegion
  path bcfCountry
  path bcfRegion
  path gatkBarcode
  path bcfBarcode
  path gatkContinent
  path bcfContinent

  output:
  path "*prediction_withMetadata.txt"
  path "*gatkcontinent_prediction_withMetadata.txt", emit: continentPredictionMerge

  script:
  """
  
  Rscript $baseDir/python_R_scripts/mergeMetada_prediction_Pfal.R ${gatkCountry} ${gatkRegion}  ${bcfCountry}  ${bcfRegion}  ${gatkBarcode}  ${bcfBarcode} ${params.metadata} ${gatkContinent} ${bcfContinent}
  """
}


process finalReport {
  // label 'rMerging'
  label 'reporting'
  publishDir "${params.reportFolder}", pattern: "*html", mode: "copy"

  input:
  path gatkCountry
  path reportableSNPs

  output:
  path "*html"

  script:
  """
  cp $baseDir/python_R_scripts/extractVOI_pfal.R .
  cp $baseDir/python_R_scripts/Pfalciparum_markdown.Rmd .
  cp $baseDir/python_R_scripts/Rplots.pdf .
  Rscript extractVOI_pfal.R  Pfalciparum_markdown.Rmd  ${reportableSNPs} ${gatkCountry} Rplots.pdf
  """
}




  // Rscript $baseDir/python_R_scripts/extractVOI__pfal_state.R  $baseDir/python_R_scripts/Pfalciparum_markdown.Rmd ${gatkCountry} ${reportableSNPs} $baseDir/python_R_scripts/Rplots.pdf