#!/bin/bash
# Run BQSR on existing markdup.bam files.
# Use this when alignment + markdup is already done but BQSR failed or was skipped.
# Expects: <OUTDIR>/<sample>/markdup.bam to exist for each sample.
set -euo pipefail

module load SAMtools/1.18-GCC-12.3.0
module load GATK/4.2.4.1-GCCcore-10.2.0-Java-1.8

GENOME="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.fa"
OUTDIR="/srv/home/aste0033/projects/MPNST/scDNA/10X/bwa_output"
JAVA_MEM="-Xmx64g"

DBSNP="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS="/mnt/iribhm/genomes/hg38-gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_INDELS="/mnt/iribhm/genomes/hg38-gatk/Homo_sapiens_assembly38.known_indels.vcf.gz"

PARALLEL=4

bqsr_sample() {
    local sample=$1
    local sampledir="${OUTDIR}/${sample}"
    local markdup_bam="${sampledir}/markdup.bam"
    local recal_table="${sampledir}/recal.table"
    local final_bam="${sampledir}/possorted_bam.bam"

    if [[ ! -f "$markdup_bam" ]]; then
        echo "[${sample}] ERROR: markdup.bam not found, skipping"
        return 1
    fi
    if [[ -f "$final_bam" ]]; then
        echo "[${sample}] possorted_bam.bam already exists, skipping"
        return 0
    fi

    echo "[${sample}] Running BaseRecalibrator..."
    gatk --java-options "$JAVA_MEM" BaseRecalibrator \
        -I "$markdup_bam" \
        -R "$GENOME" \
        --known-sites "$DBSNP" \
        --known-sites "$MILLS" \
        --known-sites "$KNOWN_INDELS" \
        -O "$recal_table" \
        --tmp-dir "$sampledir" \
        2>"${sampledir}/bqsr_recal.log"

    echo "[${sample}] Applying BQSR..."
    gatk --java-options "$JAVA_MEM" ApplyBQSR \
        -I "$markdup_bam" \
        -R "$GENOME" \
        --bqsr-recal-file "$recal_table" \
        -O "$final_bam" \
        --tmp-dir "$sampledir" \
        2>"${sampledir}/bqsr_apply.log"

    samtools index "$final_bam"
    rm "$markdup_bam" "${markdup_bam}.bai"
    echo "[${sample}] Done → ${final_bam}"
}

_pids=()
_submit() {
    bqsr_sample "$@" &
    _pids+=($!)
    if [[ ${#_pids[@]} -ge $PARALLEL ]]; then
        for _pid in "${_pids[@]}"; do
            wait "$_pid" || echo "WARNING: job $_pid exited with error"
        done
        _pids=()
    fi
}

for sample in FIT208A3 FIT208A4 FIT208A5 FIT208A6 FIT208A7 FIT208A8; do
    _submit "$sample"
done

for _pid in "${_pids[@]}"; do
    wait "$_pid" || echo "WARNING: job $_pid exited with error"
done

echo "BQSR complete for all samples."
