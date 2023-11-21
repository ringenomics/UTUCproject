cd /temp_data/DNA_data
conda activate pipeline_macvs

for sample in $(ls -d */ | cut -f1 -d'/')
do
    if [ ! -e /temp_data/DNA_data/${sample}/${sample}.dedup.numAligned.txt ]
    then
        cd ${sample}
        genefuse -r /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa -f /home/rin/Desktop/Reference_files/hg19_files_goldenpath/cancer.hg19.csv -1 ${sample}_1.fastq -2 ${sample}_2.fastq -h ${sample}_report.html -t 8 > ${sample}_genefususion_result.txt
        genefuse -r /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa -f /home/rin/Desktop/Reference_files/hg19_files_goldenpath/cancer.hg19.csv -u 1 -1 ${sample}_1.fastq -2 ${sample}_2.fastq -h ${sample}_uset1_report.html -t 8 > ${sample}_uset1_genefususion_result.txt

        trim_galore --fastqc --paired -j 25 ${sample}_1.fastq ${sample}_2.fastq

        bwa mem -R '@RG\tID:1\tSM:'${sample}'\tPL:illumina\tLB:lib1\tPU:unit1' -t 30 /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa ${sample}_1_val_1.fq ${sample}_2_val_2.fq | samtools view -bS -F 0x900 - > ${sample}.aligned.bam
        samtools view -c -F 0x4 ${sample}.aligned.bam > ${sample}.numAligned.txt

        rm -R ${sample}_1_val_1.fq ${sample}_2_val_2.fq

        samtools sort -@ 25 -m 4G -o ${sample}.aligned.bam ${sample}.aligned.bam
        samtools index ${sample}.aligned.bam

        picard MarkDuplicates --INPUT ${sample}.aligned.bam --OUTPUT ${sample}.dedup.aligned.bam --METRICS_FILE ${sample}.metrics.txt --ASSUME_SORTED True --REMOVE_DUPLICATES True
        samtools sort -@ 25 -m 4G -o ${sample}.dedup.aligned.bam ${sample}.dedup.aligned.bam
        samtools index ${sample}.dedup.aligned.bam

        samtools view -c -F 0x4 ${sample}.dedup.aligned.bam > ${sample}.dedup.numAligned.txt
        cd ..
    fi
done


