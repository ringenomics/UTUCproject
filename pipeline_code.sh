cd /Volumes/UTUCproject/Pipeline_outputs
for sample in TL-22-JICYR8PP_N_DSQ1 TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1
do
    if [ ! -e /Volumes/UTUCproject/Pipeline_outputs/${sample}/${sample}.dedup.numAligned.txt ]
    then
        cd ${sample}

        trim_galore --fastqc --paired -j 10 ${sample}_1.fastq ${sample}_2.fastq

        bwa mem -R '@RG\tID:1\tSM:'${sample}'\tPL:illumina\tLB:lib1\tPU:unit1' -t 15 /Volumes/UTUCproject/Pipeline_outputs/hg19.fa ${sample}_1_val_1.fq ${sample}_2_val_2.fq | samtools view -bS -F 0x900 - > ${sample}.aligned.bam
        samtools view -c -F 0x4 ${sample}.aligned.bam > ${sample}.numAligned.txt

        rm -R ${sample}_1_val_1.fq ${sample}_2_val_2.fq

        mkdir tmp
        samtools sort -@ 10 -m 5G -T tmp/aln.sorted -o ${sample}.aligned.bam ${sample}.aligned.bam
        samtools index ${sample}.aligned.bam

        picard MarkDuplicates --INPUT ${sample}.aligned.bam --OUTPUT ${sample}.dedup.aligned.bam --METRICS_FILE ${sample}.metrics.txt --ASSUME_SORTED True --REMOVE_DUPLICATES True
        samtools sort -@ 10 -m 3G -T tzmp/aln.sorted -o ${sample}.dedup.aligned.bam ${sample}.dedup.aligned.bam
        samtools index ${sample}.dedup.aligned.bam

        samtools view -c -F 0x4 ${sample}.dedup.aligned.bam > ${sample}.dedup.numAligned.txt
        cd ..

        #put SVIM stuff here
    fi
done