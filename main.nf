nextflow.enable.dsl=2

/*
 * Input parameters: read pairs
 * Params are stored in the params.config file
 */

// import modules for P falciparum mars

include { Trim_reads; } from './modules/trim_read'
include { PreFastqC; pre_multiqc; PostFastqC; multiqc} from './modules/fastqc'
include {BWA_index} from './modules/index'
include {BWA_align} from './modules/align'
include {Sam_sort; Picard_add_read; VCF_call; Get_Bed} from './modules/vcf_call'
include {buildsnpeff_db; annotation; vartype  } from './modules/annotation'
include {vcf_to_DF; csv_merge} from './modules/csv_merge'
include {getcoverage; WT_cov ; Trim_Stats; Reads_merge} from './modules/coverage'
include {Snpfilter; Summary_merge; Summary; Dataviz_Reportable_snps; DataViz_Novel_snps; Introns_merge } from './modules/final_snp'

//modules for WGS
include { trimReads; humanMap; noHumanReads; map3D7; mapToRefs_bwa; mapToRefs_bowtie2; fastqc_afterTrim; fastqc_beforeTrim; getDepth } from './modules/mapping.nf'
include { trimmomaticTrim; humanMap_bwa; processBAM; makeGVCFs; extractMappedBAM; genomicsDBI } from './modules/gatk4_VCF.nf'

//Import modules for geo classifier
include { trainBALK; mapPhix; fastqc_beforeTrim; fastqc_afterTrim; trimReads; humanMap; noHumanReads; indexRefs; mapToRefs; mapToRefs_bowtie2; cleanSAM_bowtie2; gatkHapCaller_bowtie2; indexBAM_bowtie2; cleanSAM; markDupes; gatkHapCaller; bcftoolsMpileup; bcftoolsCall; bcftoolsFilter_Index; gatherGVCFs; genotypeGATK; variantsToTable_gatk; gatkBarcode; predictGeo_gatk;  mergeBCFtools; variantsToTable_bcftools; bcftoolsBarcode; predictGeo_bcftools; mergePredictions; indexBAM; finalReport  } from './modules/pfal_geo_modules'




workflow P_falciparum_MaRS {

    myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")

    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    // read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).filter{ it.size()>0 }
    adapter_ch = Channel.fromPath(params.adapters, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    gff_ch = Channel.fromPath( params.gff, checkIfExists: true )
    voi_ch = Channel.fromPath( params.voi, checkIfExists: true )
  //  bed_ch = Channel.fromPath( params.bed, checkIfExists: true )


    buildsnpeff_db()
    PreFastqC(myReads)
    Trim_reads(adapter_ch.combine(myReads))
    PostFastqC(Trim_reads.out.Trimmed_fastq)
    BWA_index(ref_ch)
    BWA_align(BWA_index.out.collect().combine(Trim_reads.out.Trimmed_fastq))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    VCF_call(BWA_index.out.combine(bed_ch).combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants), buildsnpeff_db.out.buildDB)
    vartype(annotation.out.var_annotation)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    Trim_Stats(Trim_reads.out.Trimmed_stats.join(getcoverage.out.Samtools_coverage))
    WT_cov(ref_ch.combine(bed_ch).combine(voi_ch).combine(getcoverage.out.samtools_depth))
    Reads_merge(Trim_Stats.out.Reads_cov.collect())
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    Snpfilter(voi_ch.combine(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage)))
    Summary_merge(Snpfilter.out.snp_report.collect())
    Summary(Summary_merge.out)
    Introns_merge(csv_merge.out.Introns_file.collect())
    Dataviz_Reportable_snps(Summary.out.Reportable_snps)
    DataViz_Novel_snps(Summary.out.Novel_snps)
  }


workflow P_falciparum_GeoClassifier_complete_pipeline {
  // myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  
  
  // fastqc_beforeTrim(myReads)
  trimReads(myReads)
  // fastqc_afterTrim(trimReads.out.trimmedFASTQ)

  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)

  mapToRefs_bowtie2(noHumanReads.out.humanRemoveFASTQ)
  // mapToRefs_bowtie2(myReads)

  cleanSAM_bowtie2(mapToRefs_bowtie2.out.bowtieMapped)
  indexBAM_bowtie2(cleanSAM_bowtie2.out.cleanedSAM)
  gatkHapCaller_bowtie2(indexBAM_bowtie2.out.indexOut)

  gatk_allFiles = gatkHapCaller_bowtie2.out.gatkOut.collect()
  gatherGVCFs(gatk_allFiles.collect())
  genotypeGATK(gatherGVCFs.out.combinedGVCFs)
  variantsToTable_gatk(genotypeGATK.out.gatkVCF)
  gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())
  predictGeo_gatk(gatkBarcode.out.barcodeOut)


  bcftoolsMpileup(indexBAM_bowtie2.out.indexOut)
  bcftoolsCall(bcftoolsMpileup.out.pileupOut)
  bcftoolsFilter_Index(bcftoolsCall.out.callOut)
  bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
  mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())
  variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)
  bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)
  predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)
  mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut, predictGeo_gatk.out.gatkContinent, predictGeo_bcftools.out.bcfContinent)

}


workflow P_falciparum_MaRS_GeoPrediction  {

  //Set up channels
      myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
      

    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    // read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).filter{ it.size()>0 }
    adapter_ch = Channel.fromPath(params.adapters, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    gff_ch = Channel.fromPath( params.gff, checkIfExists: true )
    voi_ch = Channel.fromPath( params.voi, checkIfExists: true )
  //  bed_ch = Channel.fromPath( params.bed, checkIfExists: true )

  //build snpeff db for mars
    buildsnpeff_db()

  //index the mars references
    BWA_index(ref_ch)

  //qc of fastq reads
     PreFastqC(myReads)
    Trim_reads(adapter_ch.combine(myReads))
    humanMap(Trim_reads.out.Trimmed_fastq)
    noHumanReads(humanMap.out.sortedNoHumanBAM)
    PostFastqC(noHumanReads.out.humanRemoveFASTQ)
    multiqc(PostFastqC.out.collect())


    
    // PreFastqC(read_ch)
    // Trim_reads(adapter_ch.combine(myReads))
    // PostFastqC(Trim_reads.out.Trimmed_fastq)
   
   //run mars workflow
    BWA_align(BWA_index.out.collect().combine(noHumanReads.out.humanRemoveFASTQ))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    Get_Bed (gff_ch)
    VCF_call(BWA_index.out.combine(Get_Bed.out).combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants), buildsnpeff_db.out.buildDB)
    vartype(annotation.out.var_annotation)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    Trim_Stats(Trim_reads.out.Trimmed_stats.join(getcoverage.out.Samtools_coverage))

    WT_cov(ref_ch.combine(gff_ch).combine(voi_ch).combine(getcoverage.out.samtools_depth))
    Reads_merge(Trim_Stats.out.Reads_cov.collect())
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    Snpfilter(voi_ch.combine(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage)))
    Summary_merge(Snpfilter.out.snp_report.collect())
    Summary(Summary_merge.out)
    Introns_merge(csv_merge.out.Introns_file.collect())
    Dataviz_Reportable_snps(Summary.out.Reportable_snps)
    DataViz_Novel_snps(Summary.out.Novel_snps)

  //run geo prediction workflow
    // humanMap(trimReads.out.trimmedFASTQ)

    // noHumanReads(humanMap.out.sortedNoHumanBAM)

    mapToRefs_bowtie2(noHumanReads.out.humanRemoveFASTQ)
    cleanSAM_bowtie2(mapToRefs_bowtie2.out.bowtieMapped)
    indexBAM_bowtie2(cleanSAM_bowtie2.out.cleanedSAM)
    gatkHapCaller_bowtie2(indexBAM_bowtie2.out.indexOut)

    gatk_allFiles = gatkHapCaller_bowtie2.out.gatkOut.collect()
    gatherGVCFs(gatk_allFiles.collect())
    genotypeGATK(gatherGVCFs.out.combinedGVCFs)
    variantsToTable_gatk(genotypeGATK.out.gatkVCF)
    gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())
    predictGeo_gatk(gatkBarcode.out.barcodeOut)


    bcftoolsMpileup(indexBAM_bowtie2.out.indexOut)
    bcftoolsCall(bcftoolsMpileup.out.pileupOut)
    bcftoolsFilter_Index(bcftoolsCall.out.callOut)
    bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
    mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())
    variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)
    bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)
    predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)
    mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut, predictGeo_gatk.out.gatkContinent, predictGeo_bcftools.out.bcfContinent)

  //merge outputs into reports
    finalReport(mergePredictions.out.continentPredictionMerge, Summary.out.Summary_Snp)

}


workflow P_falciparum_GeoClassifier_useCleanReads {

  humanRemovedReads = Channel.fromFilePairs("${params.fastqDir}/*.host_removed_{R1,R2}.fastq.gz")
  // humanRemovedReads = Channel.fromFilePairs("${params.cleanReads}/*.host_removed_{R1,R2}.fastq.gz")

  mapToRefs_bowtie2(humanRemovedReads)
  // mapToRefs_bowtie2(myReads)

  cleanSAM_bowtie2(mapToRefs_bowtie2.out.bowtieMapped)
  indexBAM_bowtie2(cleanSAM_bowtie2.out.cleanedSAM)
  gatkHapCaller_bowtie2(indexBAM_bowtie2.out.indexOut)

  gatk_allFiles = gatkHapCaller_bowtie2.out.gatkOut.collect()
  gatherGVCFs(gatk_allFiles.collect())
  genotypeGATK(gatherGVCFs.out.combinedGVCFs)
  variantsToTable_gatk(genotypeGATK.out.gatkVCF)
  gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())
  predictGeo_gatk(gatkBarcode.out.barcodeOut)


  bcftoolsMpileup(indexBAM_bowtie2.out.indexOut)
  bcftoolsCall(bcftoolsMpileup.out.pileupOut)
  bcftoolsFilter_Index(bcftoolsCall.out.callOut)
  bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
  mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())
  variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)
  bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)
  predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)
  mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut, predictGeo_gatk.out.gatkContinent, predictGeo_bcftools.out.bcfContinent)

/*  
  mapToRefs(humanRemovedReads)

  cleanSAM(mapToRefs.out.mappedSAM)

  markDupes(cleanSAM.out.cleanedSAM)

  indexBAM(markDupes.out.singleFile)

  gatkHapCaller(indexBAM.out.indexOut)

  indexBAM.out.indexOut.view()

  gatk_allFiles = gatkHapCaller.out.gatkOut.collect()

  gatherGVCFs(gatk_allFiles)
  genotypeGATK(gatherGVCFs.out.combinedGVCFs)
  variantsToTable_gatk(genotypeGATK.out.gatkVCF)
  gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles)

  predictGeo_gatk(gatkBarcode.out.barcodeOut)

  bcftoolsMpileup(indexBAM.out.indexOut)
  bcftoolsCall(bcftoolsMpileup.out.pileupOut)
  bcftoolsFilter_Index(bcftoolsCall.out.callOut)


  bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()

  mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())

  variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)

  bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)

  predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)

  mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut, predictGeo_gatk.out.gatkContinent, predictGeo_bcftools.out.bcfContinent)
*/
}

workflow P_falciparum_GeoClassifier_useExistingVCF {

  gatk_allFiles = Channel.fromPath("${params.variantFolder}/*gatk_bowtie2.g.vcf")

  gatherGVCFs(gatk_allFiles.collect())

  genotypeGATK(gatherGVCFs.out.combinedGVCFs)
  variantsToTable_gatk(genotypeGATK.out.gatkVCF)
  gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())

  predictGeo_gatk(gatkBarcode.out.barcodeOut)

  bcfFiles =  Channel.fromPath("${params.variantFolder}/*bcf.vcf")

  bcftoolsFilter_Index(bcfFiles)

  bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
 
  mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())

  variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)

  bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)

  predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)

  mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut, predictGeo_gatk.out.gatkContinent, predictGeo_bcftools.out.bcfContinent)

}


//Workflow to trim reads and remove human reads. 
workflow qcReads {
  myReads =Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")

  mapPhix(myReads)
  fastqc_beforeTrim(myReads)

  trimReads(myReads)
  fastqc_afterTrim(trimReads.out.trimmedFASTQ)
  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)
  
}

//Training the balk classifier would only need to be done once. Right now it is set up for geo Combined but can change if necessary. Would need to go into the modules code
workflow P_falciparum_GeoCombined_TrainBALK{
  trainBALK()
}



workflow mapHuman_then3D7 {
  myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  
  // fastqc_beforeTrim(myReads)
  trimReads(myReads)
  // fastqc_afterTrim(trimReads.out.trimmedFASTQ)

  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)

  map3D7(noHumanReads.out.humanRemoveFASTQ)

  getDepth(map3D7.out.mapped3D7)

}


workflow map3D7_thenHuman {
  myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  
  // fastqc_beforeTrim(myReads)
  trimReads(myReads)
  // fastqc_afterTrim(trimReads.out.trimmedFASTQ)

  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)

  map3D7(noHumanReads.out.humanRemoveFASTQ)

}

workflow mapFlanking {

   myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  
  // fastqc_beforeTrim(myReads)
  trimReads(myReads)
  // fastqc_afterTrim(trimReads.out.trimmedFASTQ)

  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)

  mapToRefs_bowtie2(noHumanReads.out.humanRemoveFASTQ)

  mapToRefs_bwa(noHumanReads.out.humanRemoveFASTQ)

}

workflow just3D7 {
  myReads = Channel.fromFilePairs("$baseDir/cleanReads/*_{R1,R2}.fastq.gz")

  map3D7(myReads)
}

workflow myDepth {
  // myReads = Channel.fromFilePairs("$baseDir/cleanReads/*_{R1,R2}.fastq.gz")
  bamFiles =  Channel.fromPath("$baseDir/bowtie2Mapping_summary/*bam")
  bamFiles.view()

  getDepth(bamFiles)
}

workflow optimizedGATK4_VCF {
  myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  // chromosomes = Channel.fromList(['Nu_CHR01','Nu_CHR02','Nu_CHR03','Nu_CHR04','Nu_CHR05','Nu_CHR06','Nu_CHR07','Nu_CHR08','Nu_CHR09','Nu_CHR10','Nu_CHR11','Nu_CHR12','Nu_CHR13','Nu_CHR14'])
  chromosomes = Channel.fromList(['1', '2','3','4','5','6','7','8','9','10','11','12','13','14'])

  trimmomaticTrim(myReads)
  humanMap_bwa(trimmomaticTrim.out.trimmedFASTQ)
  processBAM(humanMap_bwa.out.sortedNoHumanBAM)
  extractMappedBAM(processBAM.out.sortedDupBam)
  makeGVCFs(extractMappedBAM.out.sortedDupBamPf, chromosomes)

  // chr_frequency =  [ "1": 2, "2": 2, "3": 2, "4": 2,"5": 2, "6": 2,"7": 2, "8": 2,"9": 2, "10": 2,"11": 2, "12": 2,"13": 2, "14": 2]
  // makeGVCFs.out.perChromGVCF.view()
  // makeGVCFs.out.perChromGVCF_tuple.view()
  fullGVCFsByChrom = makeGVCFs.out.perChromGVCF_tuple.groupTuple(by:1)
  genomicsDBI(fullGVCFsByChrom)
  

  // fullGVCFs.map{ sample, chrom, vcf -> tuple( groupKey(chrom, chr_frequency[chrom]), vcf)}.groupTuple().view()
  // fullGVCFs.groupTuple(by: 1).view()
}