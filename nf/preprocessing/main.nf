#!/usr/bin/env nextflow

Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
    .into { raw_fastq_to_trim; raw_fastq_to_qc }

mode = params.publish_dir_mode
reference = defineReference()

process fastqc{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode"  
  
  input:
    set sampleID, file(raw_fastq) from raw_fastq_to_trim

  output:
    file("**")
    file "qc/fastq/*_fastqc.{zip, html}" into fastqc_to_multiqc

  script:
  """
    mkdir -p qc/fastq
    fastqc $raw_fastq -o qc/fastq -t 8
  """
}

process bwa_index_ref {
  publishDir "$params.outdir/reference/bwa_GRCh38", mode:"$mode"	
  
  input:
  file fasta from file(reference.fasta)
  
  output:
    file("**")

  when: !params.skip_reference_indexing

  script:
  """
  bwa index $fasta
  """
}

process cutadapt {
  publishDir "$params.outdir/fastq/trimmed", mode:"$mode"
  
  tag "$sampleID"

  input:
    set sampleID, file(raw_fastq) from raw_fastq_to_qc

  output:
    set (
      sampleID,
      file("${sampleID}_trimmed_*.fastq.gz")
    ) into trimmed_fastq_to_align
    file "${sampleID}.trim.log" into cutadapt_to_multiqc  

  script:
  """
    cutadapt --cores 4 \
        -o ${sampleID}_trimmed_R1.fastq.gz \
        -p ${sampleID}_trimmed_R2.fastq.gz \
        ${raw_fastq} > ${sampleID}.trim.log 
  """
}

process bwa_align {
  
  tag "$sampleID"
  
  input:
    set (
      sampleID,
      file(fastq)
    ) from trimmed_fastq_to_align
    file bwa_ref from file(reference.bwa)

  output:
  set (
    sampleID,
    file("${sampleID}_srtd.bam")
  ) into bam

  script:
  """
  bwa mem ${params.bwa_cmd_args} ${bwa_ref}/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${fastq} > ${sampleID}.sam
  samtools sort -o ${sampleID}_srtd.bam ${sampleID}.sam
  """
}

process picard_MarkDuplicates{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode" 
  
  input:
    set sampleID, file(bam) from bam

  output:
    set (
      sampleID,
      file("bam/${sampleID}/${sampleID}.srtd.markdup.bam"),
    ) into (final_bam_index, final_bam_qc, bam_to_hmmcopy)
    file("qc/markDuplicates/${sampleID}.markDuplicates")

  script:
  """
    mkdir -p qc/markDuplicates
    mkdir -p bam/${sampleID}/
    java -jar \${PICARD_JAR} MarkDuplicates \
                PG=null \
                I=${bam} \
                O=bam/${sampleID}/${sampleID}.srtd.markdup.bam \
                M=qc/markDuplicates/${sampleID}.markDuplicates
  """
}

process samtools_bam_index{
  
  tag "$sampleID"
  publishDir "$params.outdir/bam/${sampleID}", mode:"$mode" 
  
  input:
    set sampleID, file(bam) from final_bam_index

  output:
    set (
      sampleID,
      file("${bam}.bai")
    ) into (final_bam_bai, bai_to_hmmcopy)

  script:
  """
    samtools index -@ 8 ${bam} ${bam}.bai
  """
}

process picard_CollectWgsMetrics{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode" 
  
  input:
    set (
      sampleID, 
      file(bam),
      file(bai)
    ) from final_bam_qc 
      .combine (final_bam_bai, by: 0)
    file fasta from file(reference.fasta)

  output:
    file("qc/CollectWgsMetrics/${sampleID}.CollectWgsMetrics") into picard_to_multiqc

  script:
  """
    mkdir -p qc/CollectWgsMetrics
    java -jar \${PICARD_JAR} CollectWgsMetrics \
      I=${bam} \
      O=qc/CollectWgsMetrics/${sampleID}.CollectWgsMetrics \
      R=${fasta} 
  """
}

process multiqc {
    
  publishDir "$params.outdir/qc", mode:"$mode"  
    
  input:
    file ('fastqc/*') from fastqc_to_multiqc.collect().ifEmpty([])
    file ('CollectWgsMetrics/*') from picard_to_multiqc.collect().ifEmpty([])
    file ('cutadapt/*') from cutadapt_to_multiqc.collect().ifEmpty([])

  output:
    file ("**")

  script:
  """
    multiqc .
  """
}

process hmm_readCounter{
  tag "$sampleID"
  publishDir "$params.outdir/wig", mode:"$mode"  

  input:
    set (
      sampleID, 
      file(bam),
      file(bai)
    ) from bam_to_hmmcopy
      .combine (bai_to_hmmcopy, by: 0)

  output:
    set (
      sampleID,
      file("${sampleID}.wig")
    ) into wig_to_ichorcna

  script:
  """
    readCounter --window 1000000 --quality 20 \
      --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
      $bam > ${sampleID}.wig
  """

}


process ichorcna {
  tag "$sampleID"
  publishDir "$params.outdir/ichorcna", mode:"$mode"  

  input:
    set (
      sampleID, 
      file(wig)
    ) from wig_to_ichorcna
  
  output:
    file("**") 

  script:
  """
    mkdir -p ${sampleID}
    Rscript \${ICHORCNA_DIR}/scripts/runIchorCNA.R \
      --id ${sampleID} \
      --WIG ${wig} \
      --ploidy "c(2,3)" \
      --normal "c(0.5,0.6,0.7,0.8,0.9)" \
      --maxCN 5 \
      --gcWig \${ICHORCNA_DIR}/inst/extdata/gc_hg38_1000kb.wig \
      --mapWig \${ICHORCNA_DIR}/inst/extdata/map_hg38_1000kb.wig \
      --centromere \${ICHORCNA_DIR}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
      --normalPanel \${ICHORCNA_DIR}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
      --includeHOMD False \
      --chrs "c(1:22, \'X\')" \
      --chrTrain "c(1:22)" \
      --estimateNormal True \
      --estimatePloidy True \
      --estimateScPrevalence True \
      --scStates "c(1,3)" \
      --txnE 0.9999 \
      --txnStrength 10000 \
      --outDir ${sampleID}
  """

}



/*
________________________________________________________________________________
                            F U N C T I O N S
________________________________________________________________________________
*/

def checkParamReturnFileReferences(item) {
    params."${item}" = params.references."${item}"
    return file(params."${item}")
}

def defineReference() {
    if (params.references.size() != 3) exit 1, """
    ERROR: Not all References needed found in configuration
    """
    return [
        'bwa'     : checkParamReturnFileReferences("bwa"),
        'fasta'    : checkParamReturnFileReferences("fasta"),
        'gtf'      : checkParamReturnFileReferences("gtf"),
    ]
}

