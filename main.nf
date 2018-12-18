#!/usr/bin/env nextflow

VERSION="0.2"

log.info "===================================================================="
log.info "GATK4 Best Practice Nextflow Pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz"
  log.info " "
  log.info "Mandatory arguments:"
  log.info "    --fastq1        FILE               Fastq(.gz) file for read1"
  log.info "    --fastq2        FILE               Fastq(.gz) file for read2"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --outdir        DIR                Output directory(default: ./Results)"
  log.info "    --samplename    STRING             Sample name(dafault: fastq1 basename)"
  log.info "    --rg            STRING             Read group tag(dafault: fastq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}


// Validate inputs
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; fasta_genotypegvcfs; fasta_variantrecalibrator_snps; fasta_variantrecalibrator_tranches }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
           .into { fai_bwa; fai_baserecalibrator; fai_haplotypecaller; fai_genotypegvcfs; fai_variantrecalibrator_snps; fai_variantrecalibrator_tranches }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_bwa; dict_baserecalibrator; dict_haplotypecaller; dict_genotypegvcfs; dict_variantrecalibrator_snps; dict_variantrecalibrator_tranches }
}
params.dbsnp = params.genome ? params.genomes[ params.genome ].dbsnp ?: false : false
if (params.dbsnp) {
    Channel.fromPath(params.dbsnp)
           .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp}" }
           .into { dbsnp; dbsnp_variantrecalibrator_snps; dbsnp_variantrecalibrator_indels }
}
params.dbsnp_idx = params.genome ? params.genomes[ params.genome ].dbsnp_idx ?: false : false
if (params.dbsnp_idx) {
    Channel.fromPath(params.dbsnp_idx)
           .ifEmpty { exit 1, "dbsnp_idx annotation file not found: ${params.dbsnp_idx}" }
           .into { dbsnp_idx; dbsnp_idx_variantrecalibrator_snps; dbsnp_idx_variantrecalibrator_indels }
}
params.golden_indel = params.genome ? params.genomes[ params.genome ].golden_indel ?: false : false
if (params.golden_indel) {
    Channel.fromPath(params.golden_indel)
           .ifEmpty { exit 1, "golden_indel annotation file not found: ${params.golden_indel}" }
           .into { golden_indel; golden_indel_variantrecalibrator_indels }
}
params.golden_indel_idx = params.genome ? params.genomes[ params.genome ].golden_indel_idx ?: false : false
if (params.golden_indel_idx) {
    Channel.fromPath(params.golden_indel_idx)
           .ifEmpty { exit 1, "golden_indel_idx annotation file not found: ${params.golden_indel_idx}" }
           .into { golden_indel_idx; golden_indel_idx_variantrecalibrator_indels}
}
params.hapmap_gz = params.genome ? params.genomes[ params.genome ].hapmap_gz ?: false : false
if (params.hapmap_gz) {
    Channel.fromPath(params.hapmap_gz)
           .ifEmpty { exit 1, "hapmap_gz annotation file not found: ${params.hapmap_gz}" }
           .set { hapmap_gz }
}
params.hapmap_idx_gz = params.genome ? params.genomes[ params.genome ].hapmap_idx_gz ?: false : false
if (params.hapmap_idx_gz) {
    Channel.fromPath(params.hapmap_idx_gz)
           .ifEmpty { exit 1, "hapmap_idx_gz annotation file not found: ${params.hapmap_idx_gz}" }
           .set { hapmap_idx_gz }
}
params.omni_gz = params.genome ? params.genomes[ params.genome ].omni_gz ?: false : false
if (params.omni_gz) {
    Channel.fromPath(params.omni_gz)
           .ifEmpty { exit 1, "omni_gz annotation file not found: ${params.omni_gz}" }
           .set { omni_gz }
}
params.omni_idx_gz = params.genome ? params.genomes[ params.genome ].omni_idx_gz ?: false : false
if (params.omni_idx_gz) {
    Channel.fromPath(params.omni_idx_gz)
           .ifEmpty { exit 1, "omni_idx_gz annotation file not found: ${params.omni_idx_gz}" }
           .set { omni_idx_gz }
}
params.phase1_snps = params.genome ? params.genomes[ params.genome ].phase1_snps ?: false : false
if (params.phase1_snps) {
    Channel.fromPath(params.phase1_snps)
           .ifEmpty { exit 1, "phase1_snps annotation file not found: ${params.phase1_snps}" }
           .set { phase1_snps }
}
params.phase1_snps_idx = params.genome ? params.genomes[ params.genome ].phase1_snps_idx ?: false : false
if (params.phase1_snps_idx) {
    Channel.fromPath(params.phase1_snps_idx)
           .ifEmpty { exit 1, "phase1_snps_idx annotation file not found: ${params.phase1_snps_idx}" }
           .set { phase1_snps_idx }
}
params.bwa_index_amb = params.genome ? params.genomes[ params.genome ].bwa_index_amb ?: false : false
if (params.bwa_index_amb) {
    Channel.fromPath(params.bwa_index_amb)
           .ifEmpty { exit 1, "bwa_index_amb annotation file not found: ${params.bwa_index_amb}" }
           .set { bwa_index_amb }
}
params.bwa_index_ann = params.genome ? params.genomes[ params.genome ].bwa_index_ann ?: false : false
if (params.bwa_index_ann) {
    Channel.fromPath(params.bwa_index_ann)
           .ifEmpty { exit 1, "bwa_index_ann annotation file not found: ${params.bwa_index_ann}" }
           .set { bwa_index_ann }
}
params.bwa_index_bwt = params.genome ? params.genomes[ params.genome ].bwa_index_bwt ?: false : false
if (params.bwa_index_bwt) {
    Channel.fromPath(params.bwa_index_bwt)
           .ifEmpty { exit 1, "bwa_index_bwt annotation file not found: ${params.bwa_index_bwt}" }
           .set { bwa_index_bwt }
}
params.bwa_index_pac = params.genome ? params.genomes[ params.genome ].bwa_index_pac ?: false : false
if (params.bwa_index_pac) {
    Channel.fromPath(params.bwa_index_pac)
           .ifEmpty { exit 1, "bwa_index_pac annotation file not found: ${params.bwa_index_pac}" }
           .set { bwa_index_pac }
}
params.bwa_index_sa = params.genome ? params.genomes[ params.genome ].bwa_index_sa ?: false : false
if (params.bwa_index_sa) {
    Channel.fromPath(params.bwa_index_sa)
           .ifEmpty { exit 1, "bwa_index_sa annotation file not found: ${params.bwa_index_sa}" }
           .set { bwa_index_sa }
}



/*
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */
if(params.singleEnd){
   reads="${params.reads_folder}/*.${params.reads_extension}"
  Channel
      .fromPath(reads)
      .map { file -> tuple(file.baseName, file) }
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .combine(fasta_bwa)
      .dump(tag:'input')
      .into { reads_samplename; reads_bwa }

} else if (params.pairedEnd){
  reads="${params.reads_folder}/${params.reads_prefix}_{1,2}.${params.reads_extension}"
  Channel
      .fromFilePairs(reads, size: 2)
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .combine(fastaChannel)
      .dump(tag:'input')
      .into { reads_samplename; reads_bwa }
} else {
  exit 1, "Please specify either --singleEnd or --pairedEnd to execute the pipeline!"
}

reads_samplename.subscribe { println it }



process gunzip_hapmap {
  tag "$hapmap_gz"
	publishDir "${params.outdir}/reference"

  input:
  file hapmap_gz from hapmap_gz
  file hapmap_idx_gz from hapmap_idx_gz

	output:
	file "*.vcf" into hapmap
	file "*.vcf.idx" into hapmap_idx

	"""
	gunzip -d --force $hapmap_gz
	gunzip -d --force $hapmap_idx_gz
	"""
}

process gunzip_omni {
  tag "$omni_gz"
	publishDir "${params.outdir}/reference"

  input:
  file omni_gz from omni_gz
  file omni_idx_gz from omni_idx_gz

	output:
	file "*.vcf" into omni
	file "*.vcf.idx" into omni_idx

	"""
	gunzip -d --force $omni_gz
	gunzip -d --force $omni_idx_gz
	"""
}

bwa_index = bwa_index_amb.merge(bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
bwa = reads_bwa.combine(bwa_index)

process BWA {
  tag "$reads"
	publishDir "${params.outdir}/MappedRead"
	container 'kathrinklee/bwa:latest'

	input:
  set val(name), file(reads), file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa) from bwa

	output:
	set val(name), file("${name}.sam") into sam

	"""
	bwa mem -M -R '@RG\\tID:${name}\\tSM:${name}\\tPL:Illumina' $fasta $reads > ${name}.sam
	"""

}

process BWA_sort {
  tag "$sam"
	publishDir "${params.outdir}/MappedRead"
	container 'comics/samtools:latest'

	input:
  set val(name), file(sam) from sam

	output:
	set val(name), file("${name}-sorted.bam") into bam_sort

	"""
	samtools sort -o ${name}-sorted.bam -O BAM $sam
	"""

}

process MarkDuplicates {
  tag "$bam_sort"
	publishDir "${params.outdir}/MappedRead"
	container 'broadinstitute/gatk'

	input:
	set val(name), file(bam_sort) from bam_sort

	output:
	set val(name), file("${name}_MarkDup.bam") into bam_markdup_baserecalibrator, bam_markdup_applybqsr

	"""
	gatk MarkDuplicates -I $bam_sort -M metrics.txt -O ${name}_MarkDup.bam
	"""

}

baserecalibrator_index = fasta_baserecalibrator.merge(fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx)
baserecalibrator = bam_markdup_baserecalibrator.combine(baserecalibrator_index)

process BaseRecalibrator {
  tag "$bam_markdup"
	publishDir "${params.outdir}/BaseRecalibrator"
	container 'broadinstitute/gatk:latest'

	input:
  set val(name), file(bam_markdup), file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) from baserecalibrator

	output:
	set val(name), file("${name}_recal_data.table") into baserecalibrator_table

	"""
	gatk BaseRecalibrator \
	-I $bam_markdup \
	--known-sites $dbsnp \
	--known-sites $golden_indel \
	-O ${name}_recal_data.table \
	-R $fasta
	"""
}

applybqsr = baserecalibrator_table.join(bam_markdup_applybqsr)

process ApplyBQSR {
  tag "$baserecalibrator_table"
	publishDir "${params.outdir}/BaseRecalibrator"
	container 'broadinstitute/gatk:latest'

	input:
  set val(name), file(baserecalibrator_table), file(bam_markdup) from applybqsr

	output:
	set val(name), file("${name}_bqsr.bam") into bam_bqsr

	script:
	"""
	gatk ApplyBQSR -I $bam_markdup -bqsr $baserecalibrator_table -O ${name}_bqsr.bam
	"""
}

process IndexBam {
  tag "$bam"
  container 'lifebitai/samtools:latest'

  input:
  set val(name), file(bam) from bam_bqsr

  output:
  set val(name), file("${name}_bqsr.bam"), file("${name}_bqsr.bam.bai") into indexed_bam_bqsr

  script:
  """
  samtools index $bam
  """
}

haplotypecaller_index = fasta_haplotypecaller.merge(fai_haplotypecaller, dict_haplotypecaller)
haplotypecaller = indexed_bam_bqsr.combine(haplotypecaller_index)

process HaplotypeCaller {
  tag "$bam_bqsr"
	publishDir "${params.outdir}/HaplotypeCaller"
	container 'broadinstitute/gatk:latest'

	input:
  set val(name), file(bam_bqsr), file(bai), file(fasta), file(fai), file(dict) from haplotypecaller

	output:
	set file("${name}.g.vcf"), file("${name}.g.vcf.tbi") into haplotypecaller_gvcf

	script:
	"""
  gatk HaplotypeCaller \
  -R $fasta \
  -I $bam_bqsr \
  -O ${name}.g.vcf \
  -ERC GVCF
	"""
}

process MergeVCFs {
  tag "$haplotypecaller_gvcf"
	publishDir "${params.outdir}/MergeVCFs"
	container 'broadinstitute/gatk:latest'

	input:
  set file(haplotypecaller_gvcf), file(index) from haplotypecaller_gvcf.collect()

	output:
	set file("merged.g.vcf"), file("merged.g.vcf.tbi") into mergevcfs

	script:
	"""
  gatk MergeVcfs \
  --INPUT= $haplotypecaller_gvcf
  --OUTPUT= merged.g.vcf
	"""
}


process HardFilterVcf {
  tag "$merged_vcf"
	publishDir "${params.outdir}"
	container 'broadinstitute/gatk:latest'

	input:
  set file(merged_vcf), file(index) from mergevcfs

	output:
  set file("output.vcf"), file("output.vcf.tbi") into results

	script:
	"""
  gatk VariantFiltration \
  -V $merged_vcf \
  --filterExpression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
  --filterName "HardFiltered" \
  -O output.vcf
	"""
}
