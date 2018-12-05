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
           .into { fasta; fasta_bwa }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
           .into { fai; fai_bwa }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict; dict_bwa }
}
params.dbsnp = params.genome ? params.genomes[ params.genome ].dbsnp ?: false : false
if (params.dbsnp) {
    Channel.fromPath(params.dbsnp)
           .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp}" }
           .set { dbsnp }
}
params.dbsnp_idx = params.genome ? params.genomes[ params.genome ].dbsnp_idx ?: false : false
if (params.dbsnp_idx) {
    Channel.fromPath(params.dbsnp_idx)
           .ifEmpty { exit 1, "dbsnp_idx annotation file not found: ${params.dbsnp_idx}" }
           .set { dbsnp_idx }
}
params.golden_indel = params.genome ? params.genomes[ params.genome ].golden_indel ?: false : false
if (params.golden_indel) {
    Channel.fromPath(params.golden_indel)
           .ifEmpty { exit 1, "golden_indel annotation file not found: ${params.golden_indel}" }
           .set { golden_indel }
}
params.golden_indel_idx = params.genome ? params.genomes[ params.genome ].golden_indel_idx ?: false : false
if (params.golden_indel_idx) {
    Channel.fromPath(params.golden_indel_idx)
           .ifEmpty { exit 1, "golden_indel_idx annotation file not found: ${params.golden_indel_idx}" }
           .set { golden_indel_idx }
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


bwa_index = Channel.from(bwa_index_amb).merge(bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa).collect()


/*
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */
if(params.singleEnd){
   reads="${params.reads_folder}/*.${params.reads_extension}"
  Channel
      .fromFilePairs( reads, size: params.singleEnd ? 1 : 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .dump(tag:'input')
      .into { reads_samplename; reads_bwa }

} else if (params.pairedEnd){
  reads="${params.reads_folder}/${params.reads_prefix}_{1,2}.${params.reads_extension}"
  Channel
      .fromFilePairs(reads, size: 2)
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .dump(tag:'input')
      .into { reads_samplename; reads_bwa }
} else {
  exit 1, "Please specify either --singleEnd or --pairedEnd to execute the pipeline!"
}

reads_samplename.first().subscribe { println it }



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

process BWA {
  tag "$reads"
	publishDir "${params.outdir}/MappedRead"
	container 'kathrinklee/bwa:latest'

	input:
	file fasta from fasta_bwa
	file bwa_index from bwa_index
	file reads from reads_bwa

	output:
	file 'aln-pe.sam' into samfile

	"""
	bwa mem -M -R '@RG\\tID:${params.reads_prefix}\\tSM:${params.reads_prefix}\\tPL:Illumina' $fasta *.${params.reads_extension} > aln-pe.sam
	"""

}

// process BWA_sort {
// 	publishDir "${params.outdir}/MappedRead"
// 	container 'comics/samtools:latest'
//
// 	input:
// 	file samfile
//
// 	output:
// 	file 'aln-pe-sorted.bam' into bam_sort
//
// 	"""
// 	samtools sort -o aln-pe-sorted.bam -O BAM $samfile
// 	"""
//
// }
//
// process MarkDuplicates {
// 	publishDir "${params.outdir}/MappedRead"
// 	container 'broadinstitute/gatk'
//
// 	input:
// 	file bam_sort
//
// 	output:
// 	file 'aln-pe_MarkDup.bam' into bam_markdup
//
// 	"""
// 	gatk MarkDuplicates -I $bam_sort -M metrics.txt -O aln-pe_MarkDup.bam
// 	"""
//
// }
//
// process BaseRecalibrator {
// 	publishDir "${params.outdir}/BaseRecalibrator"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file reference
// 	file reference_fai
// 	file reference_dict
// 	file bam_markdup
// 	file dbsnp
// 	file dbsnp_idx
// 	file golden_indel
// 	file golden_indel_idx
//
// 	output:
// 	file 'recal_data.table' into BaseRecalibrator_table
//
// 	"""
// 	gatk BaseRecalibrator \
// 	-I $bam_markdup \
// 	--known-sites $dbsnp \
// 	--known-sites $golden_indel \
// 	-O recal_data.table \
// 	-R $reference
// 	"""
// }
//
// process ApplyBQSR {
// 	publishDir "${params.outdir}/BaseRecalibrator"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file BaseRecalibrator_table
// 	file bam_markdup
//
// 	output:
// 	file 'aln-pe_bqsr.bam' into bam_bqsr
//
// 	script:
// 	"""
// 	gatk ApplyBQSR -I $bam_markdup -bqsr $BaseRecalibrator_table -O aln-pe_bqsr.bam
// 	"""
// }
//
// process HaplotypeCaller {
// 	publishDir "${params.outdir}/HaplotypeCaller"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file reference
// 	file reference_fai
// 	file reference_dict
// 	file bam_bqsr
//
// 	output:
// 	file 'haplotypecaller.g.vcf' into haplotypecaller_gvcf
//
// 	script:
// 	"""
// 	gatk HaplotypeCaller -I $bam_bqsr -O haplotypecaller.g.vcf --emit-ref-confidence GVCF -R $reference
// 	"""
// }
//
// process GenotypeGVCFs {
// 	publishDir "${params.outdir}/HaplotypeCaller"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file reference
// 	file reference_fai
// 	file reference_dict
// 	file haplotypecaller_gvcf
//
// 	output:
// 	file 'haplotypecaller.vcf' into haplotypecaller_vcf
//
// 	script:
// 	"""
// 	gatk GenotypeGVCFs --variant haplotypecaller.g.vcf -R $reference -O haplotypecaller.vcf
// 	"""
// }
//
//
//
// process VariantRecalibrator_SNPs {
// 	publishDir "${params.outdir}/VariantRecalibrator"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file reference
// 	file reference_fai
// 	file reference_dict
// 	file haplotypecaller_vcf
// 	file hapmap
// 	file hapmap_idx
// 	file omni
// 	file omni_idx
// 	file phase1_snps
// 	file phase1_snps_idx
// 	file dbsnp
// 	file dbsnp_idx
//
// 	output:
// 	file 'recalibrate_SNP.recal' into variantrecalibrator_recal
// 	file 'recalibrate_SNP.recal.idx' into variantrecalibrator_recal_idx
// 	file 'recalibrate_SNP.tranches' into variantrecalibrator_tranches
//
// 	script:
// 	"""
// 	gatk VariantRecalibrator \
// 	-V $haplotypecaller_vcf \
//  	-R $reference \
// 	-resource hapmap,known=false,training=true,truth=true,prior=15.0:./$hapmap \
// 	-resource omni,known=false,training=true,truth=true,prior=12.0:./$omni \
//     	-resource 1000G,known=false,training=true,truth=false,prior=10.0:./$phase1_snps \
//     	-resource dbsnp,known=true,training=false,truth=false,prior=2.0:./$dbsnp \
// 	-an DP \
//     	-an QD \
// 	-an FS \
//     	-an SOR \
//     	-an MQ \
//     	-an MQRankSum \
//     	-an ReadPosRankSum \
//     	-mode SNP \
//     	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
// 	--max-gaussians 8 \
//     	-O recalibrate_SNP.recal \
//     	--tranches-file recalibrate_SNP.tranches \
// 	"""
// }
//
//
//
// process ApplyVQSR_SNPs {
// 	publishDir "${params.outdir}/VariantRecalibrator"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file haplotypecaller_vcf
// 	file variantrecalibrator_recal
// 	file variantrecalibrator_recal_idx
// 	file variantrecalibrator_tranches
//
// 	output:
// 	file 'recalibrated_snps_raw_indels.vcf' into recalibrated_snps_raw_indels
//
// 	script:
// 	"""
// 	gatk ApplyVQSR \
// 	-V $haplotypecaller_vcf \
// 	--recal-file $variantrecalibrator_recal \
// 	--tranches-file $variantrecalibrator_tranches \
// 	-mode SNP \
// 	-ts-filter-level 99.0 \
// 	-O recalibrated_snps_raw_indels.vcf
// 	"""
// }
//
//
//
// process VariantRecalibrator_INDELs {
// 	publishDir "${params.outdir}/VariantRecalibrator"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file reference
// 	file reference_fai
// 	file reference_dict
// 	file recalibrated_snps_raw_indels
// 	file dbsnp
// 	file dbsnp_idx
// 	file golden_indel
// 	file golden_indel_idx
//
// 	output:
// 	file 'recalibrate_INDEL.recal' into variantrecalibrator_indel_recal
// 	file 'recalibrate_INDEL.recal.idx' into variantrecalibrator_indel_recal_idx
// 	file 'recalibrate_INDEL.tranches' into variantrecalibrator_indel_tranches
//
// 	script:
// 	"""
// 	gatk VariantRecalibrator \
// 	-V $recalibrated_snps_raw_indels \
//  	-R $reference \
// 	--resource mills,known=false,training=true,truth=true,prior=12.0:./$golden_indel \
//     	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:./$dbsnp \
// 	-an QD \
//     	-an DP \
//     	-an FS \
// 	-an SOR \
//     	-an MQRankSum \
//     	-an ReadPosRankSum \
//     	-mode INDEL \
//     	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
// 	--max-gaussians 4 \
//     	-O recalibrate_INDEL.recal \
//     	--tranches-file recalibrate_INDEL.tranches \
// 	"""
// }
//
// process ApplyVQSR_INDELs {
// 	publishDir "${params.outdir}/VariantRecalibrator"
// 	container 'broadinstitute/gatk:latest'
//
// 	input:
// 	file recalibrated_snps_raw_indels
// 	file variantrecalibrator_indel_recal
// 	file variantrecalibrator_indel_recal_idx
// 	file variantrecalibrator_indel_tranches
//
// 	output:
// 	file 'recalibrated_variants.vcf' into recalibrated_variants_vcf
//
// 	script:
// 	"""
// 	gatk ApplyVQSR \
// 	-V $recalibrated_snps_raw_indels \
// 	--recal-file $variantrecalibrator_indel_recal \
// 	--tranches-file $variantrecalibrator_indel_tranches \
// 	-mode INDEL \
// 	-ts-filter-level 99.0 \
// 	-O recalibrated_variants.vcf
// 	"""
// }
//
// process copy {
// 	publishDir "${params.outdir}", mode: 'copy'
//
// 	input:
// 	file recalibrated_variants_vcf
//
// 	output:
// 	file "${params.samplename}.vcf"
//
// 	"""
// 	mv $recalibrated_variants_vcf ${params.samplename}.vcf
// 	"""
//
// }
