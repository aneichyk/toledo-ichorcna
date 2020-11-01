#!/usr/bin/env nextflow

if(params.readPaths){
    Channel
        .from(params.readPaths)
        .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
        .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        .into { raw_read_files; read_files_fastqc }
} else {
    Channel
        .fromFilePairs( params.reads, size: 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .into { raw_read_files; read_files_fastqc }
}

mode = params.publish_dir_mode
reference = defineReference()

/*
________________________________________________________________________________
                            P R O C E S S E S
________________________________________________________________________________
*/


process fastqc{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode"  
  
  input:
    set sampleID, file(raw_fastq) from read_files_fastqc

  output:
    file("**")

  script:
  """
    mkdir -p qc/fastq
    fastqc $raw_fastq -o qc/fastq
  """
}

process bwa_index_ref {
  publishDir "$params.outdir/reference", mode:"$mode"	
  
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

process trimgalore {
  publishDir "$params.outdir/fastq", mode:"$mode"
  
  tag "$sampleID"

  input:
    set sampleID, file(raw_fastq) from raw_read_files

  output:
  set (
    sampleID,
    file("${sampleID}_val_*.fq.gz")
  ) into trimmed_reads

  script:
  """
    trim_galore --paired ${raw_fastq}
    ls -la
  """

}

process bwa_align {
  
  tag "$sampleID"
  
  publishDir "$params.outdir/bam/bwa/", mode:"$mode"	
  
  input:
    set (
      sampleID,
      file(fastq)
    ) from trimmed_reads
    file bwa_ref from file(reference.bwaref)

  output:
  set (
    sampleID,
    file("${sampleID}_sorted.bam")
  ) into bam

  script:
  """
  bwa mem $params.bwa_cmd_args \\
    ${bwa_ref}/Homo_sapiens_assembly19.fasta \\
    $paired_fastq > ${sampleID}_paired.sam
  samtools sort -o ${sampleID}_paired_sorted.bam ${sampleID}_paired.sam
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

