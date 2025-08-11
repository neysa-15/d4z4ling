INPUT_PATH=/g/data/kr68/neysa/fshd_pipeline/benchmark_ubam_path.tsv

logs="/g/data/kr68/fshd/logs"
timestamp=$(date +%Y%m%d_%H%M%S)

while IFS=$'\t' read -r sample outdir ubam; do
    echo "Sample: $sample"
    echo "  OUTDIR: $outdir"
    echo "  uBAM: $ubam"

    qsub -v SAMPLE=$sample,OUTDIR=$outdir,UBAM=$ubam \
        -o "${logs}/${sample}.o${timestamp}" \
        -e "${logs}/${sample}.e${timestamp}" \
        /g/data/kr68/neysa/fshd_pipeline/nci_run.sh

done < "$INPUT_PATH"
