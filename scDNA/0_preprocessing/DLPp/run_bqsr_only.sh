#!/bin/bash
# Run BQSR on per-sample DLP+ BAMs when alignment + merge succeeded but
# BQSR failed or was skipped.
# Expects: <OUTDIR>/<chip>/<sample>.bam (produced by fastq_to_bam_dlpp.sh).
#
# Input:  <OUTDIR>/<chip>/<sample>.bam
# Output: <OUTDIR>/<chip>/<sample>_bqsr.bam

set -euo pipefail

module load SAMtools/1.12-GCC-10.2.0
module load GATK/4.2.4.1-GCCcore-10.2.0-Java-1.8

GENOME="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.fa"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/DLPp/bwa_output"
JAVA_MEM="-Xmx64g"

DBSNP="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS="/mnt/iribhm/genomes/hg38-gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_INDELS="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.known_indels.vcf.gz"

PARALLEL=2   # at most 2 concurrent BQSR jobs (each uses ~64 GB GATK heap)

# Sample → chip mapping
declare -A CHIP=(
    ["R1"]="LES4677A1"
    ["R2"]="LES4677A1"
    ["R3"]="LES4677A2"
    ["R4"]="LES4677A2"
    ["R5"]="LES4677A4"
    ["P"]="LES4677A4"
)
SAMPLES=(R1 R2 R3 R4 R5 P)

bqsr_sample() {
    local sample=$1
    local chip="${CHIP[$sample]}"
    local sampledir="${OUTDIR}/${chip}"
    local in_bam="${sampledir}/${sample}.bam"
    local recal_table="${sampledir}/${sample}.recal.table"
    local out_bam="${sampledir}/${sample}_bqsr.bam"

    if [[ ! -f "$in_bam" ]]; then
        echo "[${sample}] ERROR: ${in_bam} not found, skipping"
        return 1
    fi
    if [[ -f "$out_bam" ]]; then
        echo "[${sample}] ${out_bam} already exists, skipping"
        return 0
    fi

    echo "[${sample}] Running BaseRecalibrator..."
    gatk --java-options "$JAVA_MEM" BaseRecalibrator \
        -I "$in_bam" \
        -R "$GENOME" \
        --known-sites "$DBSNP" \
        --known-sites "$MILLS" \
        --known-sites "$KNOWN_INDELS" \
        -O "$recal_table" \
        --tmp-dir "$sampledir" \
        2>"${sampledir}/${sample}.bqsr_recal.log"

    echo "[${sample}] Applying BQSR..."
    gatk --java-options "$JAVA_MEM" ApplyBQSR \
        -I "$in_bam" \
        -R "$GENOME" \
        --bqsr-recal-file "$recal_table" \
        -O "$out_bam" \
        --tmp-dir "$sampledir" \
        2>"${sampledir}/${sample}.bqsr_apply.log"

    samtools index "$out_bam"
    echo "[${sample}] Done → ${out_bam}"
}

_pids=()
for sample in "${SAMPLES[@]}"; do
    bqsr_sample "$sample" &
    _pids+=($!)
    if [[ ${#_pids[@]} -ge $PARALLEL ]]; then
        for _pid in "${_pids[@]}"; do wait "$_pid" || echo "WARNING: job $_pid exited with error"; done
        _pids=()
    fi
done
for _pid in "${_pids[@]}"; do wait "$_pid" || echo "WARNING: job $_pid exited with error"; done

echo "BQSR complete."
