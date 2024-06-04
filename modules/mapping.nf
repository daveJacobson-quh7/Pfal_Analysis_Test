process trimReads {
  label 'bioinformaticProcessing'
  // label 'mappingSoftware'

  input:
  tuple val(sample_id), path(sample_files)


  output:
  tuple val(sample_id), path("*.bbduk.fastq.gz"), emit: trimmedFASTQ


  script:


   """

   bbduk.sh  -Xmx${params.bbmergeRAM} ktrim=r k=23 mink=11 hdist=1 edist=0 ref=${params.adapters} qtrim=rl trimq=20 minlength=50 trimbyoverlap=t minoverlap=24 qin=33 in=${sample_files[0]} in2=${sample_files[1]} out=${sample_id}.q20_R1.bbduk.fastq out2=${sample_id}.q20_R2.bbduk.fastq stats=${sample_id}.stats.txt

   gzip *.fastq
  
  """
}

process fastqc_beforeTrim {

  // label 'fastqc_container'

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
 publishDir "${params.reportingFolder}", pattern: "*testHuman", mode: "copy"
  label 'bioinformaticProcessing'

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*.bothReadsUnmapped_sorted.bam", emit: sortedNoHumanBAM
  path "*testHuman"


  script:
  """
  bowtie2 -p 4 -x ${params.humanRef} -1 ${sample_files[0]} -2 ${sample_files[1]}   -S ${sampleID}_mapped_and_unmapped.sam &> ${sampleID}_q20.testHuman

  echo 'after human mapping'

  samtools view -bS ${sampleID}_mapped_and_unmapped.sam > ${sampleID}_mapped_and_unmapped.bam

  samtools view -b -f 12 -F 256 ${sampleID}_mapped_and_unmapped.bam  > ${sampleID}_bothReadsUnmapped.bam 

  samtools sort -n -m 5G -@ 2 ${sampleID}_bothReadsUnmapped.bam -o ${sampleID}_q20.bothReadsUnmapped_sorted.bam

  """
}



process noHumanReads{
    
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

process map3D7 {
    publishDir "${params.reportingFolder}", pattern: "*test3D7", mode: "copy"
    publishDir "${params.reportingFolder}", pattern: "*mapped.bam", mode: "copy"
    label 'bioinformaticProcessing'

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*3d7_mapped.bam", emit: mapped3D7
  path "*test3D7"


  script:
  """
  bowtie2 -p 4 -x ${params.ref3D7} -1 ${sample_files[0]} -2 ${sample_files[1]}   -S ${sampleID}_3d7mapped_and_unmapped.sam &> ${sampleID}.test3D7

  echo 'after 3d7 mapping'

  samtools view -b -F 4 ${sampleID}_3d7mapped_and_unmapped.sam > ${sampleID}_3d7_mapped.bam
  """

}


process mapToRefs_bowtie2 {
    publishDir "${params.reportingFolder}", pattern: "*mapStats*", mode: "copy"
    publishDir "${params.mappingOutput}", pattern: "*bam", mode: "copy"
    label 'bioinformaticProcessing'

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*bam", emit: mappedBowtie2
  path "*mapStats*", emit: mapStat
 


  script:
  """
  bowtie2 -p 8 -x ${params.flankingRef_bowtie2} -1 ${sample_files[0]} -2 ${sample_files[1]}   -S ${sampleID}_mapped_and_unmapped.sam &> ${sampleID}.mapStats_${params.gene}

  echo 'after 3d7 mapping'

  samtools view -b -F 4 ${sampleID}_mapped_and_unmapped.sam > ${sampleID}_${params.gene}.bowtie2_mapped.bam
  """

}

process mapToRefs_bwa {

  // label 'mappingSoftware'
  label 'bioinformaticProcessing'
  publishDir "${params.mappingOutput}", pattern: "*bam", mode: "copy"
  

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*bwa_mapped.bam", emit: mappedBWA

  script:
  """
  bwa mem -t 16 -M -R  "@RG\\tID:'${sampleID}'\\tLB:'${sampleID}'\\tPL:illumina\\tSM:'${sampleID}'\\tPU:'${sampleID}'" ${params.flankingRef_bwa} ${sampleID}.host_removed_R1.fastq.gz ${sampleID}.host_removed_R2.fastq.gz > ${sampleID}.Map.sam

  samtools view -b -F 4 ${sampleID}.Map.sam > ${sampleID}_${params.gene}.bwa_mapped.bam

  """
}


process getDepth {

  // publishDir "${params.reportingFolder}", pattern: "*mapStats*", mode: "copy"
    publishDir "${params.reportingFolder}", pattern: "*bedtoolsCoverage.txt", mode: "copy"
    label 'bed'

  input:
  path inFile

  output:
  // path "*bam", emit: mappedBowtie2
  // path "*mapStats*", emit: mapStat
  path "*bedtoolsCoverage.txt"
 


  script:
  """

  echo ${inFile} > processingNow
  samtools sort ${inFile} > ${inFile.simpleName}_3d7_mapped_sorted.bam
  bedtools genomecov -ibam ${inFile.simpleName}_3d7_mapped_sorted.bam -g $baseDir/REFERENCES/Pf3d7_fullGenome_genomeFile.txt > ${inFile.simpleName}.default_bedtoolsCoverage.txt
  bedtools genomecov -ibam ${inFile.simpleName}_3d7_mapped_sorted.bam -g $baseDir/REFERENCES/Pf3d7_fullGenome_genomeFile.txt -bg > ${inFile.simpleName}.bedGraph_bedtoolsCoverage.txt


  """

}
