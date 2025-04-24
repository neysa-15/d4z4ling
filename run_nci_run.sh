logs="/g/data/kr68/fshd/logs"
timestamp=$(date +%Y%m%d_%H%M%S)
sample=${1}
outdir=${2}
ubam=${3}
qsub -v SAMPLE=${sample},OUTDIR=${outdir},UBAM=${ubam} -o "${logs}/${sample}.o${timestamp}" -e "${logs}/${sample}.e${timestamp}" ./nci_run.sh 
