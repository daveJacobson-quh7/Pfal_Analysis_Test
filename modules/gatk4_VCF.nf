process trimmomaticTrim {
  label 'optiGatk4'

  input:
  tuple val(sample_id), path(sample_files)


  output:
  tuple val(sample_id), path("*_paired.fq.gz"), emit: trimmedFASTQ


  script:
    """
    TrimmomaticPE ${sample_files[0]} ${sample_files[1]} ${sample_id}_R1_paired.fq.gz ${sample_id}_R1_unpaired.fq.gz ${sample_id}_R2_paired.fq.gz ${sample_id}_R2_unpaired.fq.gz  ILLUMINACLIP:${params.libKit}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36  SLIDINGWINDOW:5:20 -threads ${params.threads}

    """
}




process humanMap_bwa{
//  publishDir "${params.processingOut}", pattern: "*testHuman", mode: "copy"
  label 'optiGatk4'

  input:
  tuple val(sampleID), path(sample_files)

  output:
  path "*bam", emit: sortedNoHumanBAM
  // path "*testHuman"


  script:
  """
  bwa mem -t ${params.threads} -M -R "@RG\\tID:'${sampleID}'\\tLB:'${sampleID}'\\tPL:illumina\\tSM:'${sampleID}'\\tPU:'${sampleID}'" ${params.gatk4PipelineRefs}/Pf3D7_human.fa ${sample_files[0]} ${sample_files[1]} > ${sampleID}.sam

  gatk --java-options -Xmx${params.gatk4_RAM} SamFormatConverter -R ${params.gatk4PipelineRefs}/Pf3D7_human.fa -I ${sampleID}.sam -O ${sampleID}.bam

  """
}


process processBAM {
//  publishDir "${params.processingOut}", pattern: "*testHuman", mode: "copy"
  label 'optiGatk4'

  input:
  path sampleBam

  output:
  path "*sorted.dup.bam", emit: sortedDupBam
  path "*dup_metrics.tsv", emit: dupMetrics



  script:
  """
  gatk --java-options -Xmx${params.gatk4_RAM}   CleanSam -R ${params.gatk4PipelineRefs}/Pf3D7_human.fa  -I ${sampleBam} -O ${sampleBam.simpleName}.clean.bam
  gatk --java-options -Xmx${params.gatk4_RAM}   SortSam -R ${params.gatk4PipelineRefs}/Pf3D7_human.fa  -I ${sampleBam.simpleName}.clean.bam -O ${sampleBam.simpleName}.sorted.bam -SO coordinate --CREATE_INDEX true
  gatk --java-options -Xmx${params.gatk4_RAM}   MarkDuplicates -R ${params.gatk4PipelineRefs}/Pf3D7_human.fa  -I ${sampleBam.simpleName}.sorted.bam -O ${sampleBam.simpleName}.sorted.dup.bam -M ${sampleBam.simpleName}_dup_metrics.tsv
  
  """
}

process extractMappedBAM {
//  publishDir "${params.processingOut}", pattern: "*testHuman", mode: "copy"
  label 'optiGatk4'

  input:
  path cleanSampleBam

  output:
  // path "*sorted.dup.pf.ba*", emit: sortedDupBamPf
  path "*sorted.dup.hs.bam", emit: sortedDupBamHs
  // path "*dup_metrics.tsv", emit: dupMetrics
  tuple val("${cleanSampleBam.simpleName}"), path("*sorted.dup.pf.ba*"), emit: sortedDupBamPf

  script:
  """
  samtools view -b -h ${cleanSampleBam} -T ${params.gatk4PipelineRefs}/Pf3D7.fasta -L ${params.gatk4PipelineRefs}/Pf3D7_core.bed > ${cleanSampleBam.simpleName}.sorted.dup.pf.bam
  samtools view -b -h ${cleanSampleBam}  -T ${params.gatk4PipelineRefs}/Pf3D7_human.fa -L ${params.gatk4PipelineRefs}/human.bed > ${cleanSampleBam.simpleName}.sorted.dup.hs.bam
  samtools index -bc ${cleanSampleBam.simpleName}.sorted.dup.pf.bam

  """
}

// Come back to caluclating mapping statistics from the QC pipeline for optimal gatk4. They aren't neccessary for running the pipeline, from what I can tell. But nice to have for later on

process makeGVCFs{

  //  publishDir "${params.processingOut}", pattern: "*testHuman", mode: "copy"
  label 'optiGatk4'

  input:
  tuple val(sampleID), path(coreBamPf)
  each myChrom

  output:
  // // path "*sorted.dup.pf.ba*", emit: sortedDupBamPf
  // path "*sorted.dup.hs.bam", emit: sortedDupBamHs
  // // path "*dup_metrics.tsv", emit: dupMetrics
  // tuple val(${cleanSampleBam.simpleName}), path("*sorted.dup.pf.ba**"), emit: sortedDupBamPf
  path "${sampleID}.chr${myChrom}.g.vcf", emit: perChromGVCF
  tuple val("${sampleID}"), val("${myChrom}"), path("${sampleID}.chr${myChrom}.g.vcf"), path("${sampleID}.chr${myChrom}.g.vcf.idx"), emit: perChromGVCF_tuple


  script:
  """
  gatk   HaplotypeCaller -R ${params.gatk4PipelineRefs}/Pf3D7.fasta -I ${sampleID}.sorted.dup.pf.bam -ERC GVCF --native-pair-hmm-threads ${params.threads} -ploidy 2 -O  ${sampleID}.chr${myChrom}.g.vcf --assembly-region-padding 100 --max-num-haplotypes-in-population 128 --kmer-size 10 --kmer-size 25 --min-dangling-branch-length 4 --heterozygosity 0.0029 --indel-heterozygosity 0.0017 --min-assembly-region-size 100 -L ${params.gatk4PipelineRefs}/core_chr${myChrom}.list -mbq 5 -DF MappingQualityReadFilter --base-quality-score-threshold 12

  """

}
  // echo "${vcfList.join('\n')}" > vcf.list
process genomicsDBI {
 label 'optiGatk4'

  input:
  tuple val(samples), val(chromIn), path(vcfList), path(indexList)

  output:
  // tuple val("gatk"), path("gatk.g.vc*"), emit: combinedGVCFs
  tuple val("${chromIn}"), path("chr${chromIn}_database"), emit: outDB

  script:
  """
  echo ${samples} > sampleList.txt
  echo ${chromIn} > chromsAnalyzed.txt
  echo ${vcfList} > myVCFs.txt
  echo "${vcfList.join('\n')}" > vcf.list
  

  gatk --java-options -Xmx${params.gatk4_RAM} GenomicsDBImport \
    --variant vcf.list \
    --genomicsdb-workspace-path chr${chromIn}_database \
    --intervals ${params.gatk4PipelineRefs}/core_chr${chromIn}.list --batch-size 100 \
    --reader-threads 24 --genomicsdb-segment-size 8048576 \
    --genomicsdb-vcf-buffer-size 160384


  """
}

//need to find a way to do the regions in each file - 5kb is good
//Maybe just cat in a loop? Come back to splittying up?
process genotypeGVCFs {
 label 'optiGatk4'

  input:
  tuple val(samples), val(chromIn), path(vcfList), path(indexList)
  path chrDB

  // output:
  // tuple val("gatk"), path("gatk.g.vc*"), emit: combinedGVCFs

  script:
  """
  gatk --java-options -Xmx${params.gatk4_RAM}  GenotypeGVCFs --genomicsdb-use-bcf-codec true  -R ${params.gatk4PipelineRefs}/Pf3D7.fasta -V gendb://${chrDB} --max-genotype-count 1024 -O chr"$i"_part"$j".vcf.gz --tmp-dir ${params.gatk4PipelineRefs} -stand-call-conf 30 -L "$j"



  """
}
/*
  gatk CombineGVCFs -R ${params.ampliseqRef} --variant vcf.list -O gatk.g.vcf

  gatk --java-options -Xmx${params.gatk4_RAM} GenomicsDBImport \
    --sample-name-map gvcf_chr"$i"_list.tsv \
    --genomicsdb-workspace-path chr${myChrom}_database \
    --intervals ${params.gatk4PipelineRefs}/core_chr${myChrom}.list --batch-size 100 \
    --reader-threads 24 --genomicsdb-segment-size 8048576 \
    --genomicsdb-vcf-buffer-size 160384
*/
