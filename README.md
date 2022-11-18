# CHIP-seq

01. Trimming

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=trim_galore
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB

ARRAY=(WT_18DAP_embryo_H3K9me2_ChIP_rep1 WT_18DAP_embryo_H3K9me2_ChIP_rep2 WT_18DAP_embryo_H3K9me2_Input_rep1 WT_18DAP_embryo_H3K9me2_Input_rep2)

trim_galore --illumina --fastqc --paired ../../Raw_Reads/DDM1_ChIP_seq_PC_paper_raw_reads/${ARRAY[$SLURM_ARRAY_TASK_ID]}_R1.fastq.gz \
../../Raw_Reads/DDM1_ChIP_seq_PC_paper_raw_reads/${ARRAY[$SLURM_ARRAY_TASK_ID]}_R2.fastq.gz


02. Mapping 

#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=Bowtie2_ChIP_mapping
#SBATCH --cpus-per-task=16
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH -t 60:00:00

ARRAY=(WT_18DAP_embryo_H3K9me2_Input_rep1 WT_18DAP_embryo_H3K9me2_Input_rep2)

bowtie2 -p 16 --no-unal --no-discordant -x ../../Index/Bowtie2_maize_index/Zm_V5_DNA.toplevel -1 ../../Clean_Reads/DDM1_ChIP_seq_PC_paper_Clean_reads/${ARRAY[$SLURM_ARRAY_TASK_ID]}_R1_val_1.fq.gz \
-2 ../../Clean_Reads/DDM1_ChIP_seq_PC_paper_Clean_reads/${ARRAY[$SLURM_ARRAY_TASK_ID]}_R2_val_2.fq.gz |samtools sort -@ 2 -O Bam -o ${ARRAY[$SLURM_ARRAY_TASK_ID]}.sorted.bam -



03. Remove_duplicates
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=rm_PCR_duplicates
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30GB
#SBATCH -t 12:00:00


java -jar /data6/tool/picard-tools-2.10.10/picard.jar MarkDuplicates I=WT_18DAP_embryo_IgG_ChIP_rep1.sorted.bam O=WT_18DAP_embryo_IgG_ChIP_rep1.sorted_rmdup.bam \
M=WT_18DAP_embryo_IgG_ChIP_rep1.merics.txt AS=true REMOVE_DUPLICATES=true

java -jar /data6/tool/picard-tools-2.10.10/picard.jar MarkDuplicates I=WT_18DAP_embryo_IgG_ChIP_rep2.sorted.bam O=WT_18DAP_embryo_IgG_ChIP_rep2.sorted_rmdup.bam \
M=WT_18DAP_embryo_IgG_ChIP_rep2.merics.txt AS=true REMOVE_DUPLICATES=true


04. Filtering (Optional)

samtools view -@ 3 -b -q 5 ../ZmAGO4_ChIP_rep2.sorted_rmdup.bam -o ZmAGO4_ChIP_rep2.sorted_rmdup_MAPQ5.bam


05. Calling peaks by macs2
#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=macs2_callpeaks
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1

macs2 callpeak -B --SPMR -t ZmDDM1A_ChIP_rep1.sorted_rmdup.bam -c ZmDDM1A_Input_rep1.sorted_rmdup.bam -q 0.01 --outdir ZmDDM1A_rep1_ChIP_output


