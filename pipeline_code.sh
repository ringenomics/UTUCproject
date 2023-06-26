#Run in pipeline env

cd /Volumes/UTUCproject/DNA_Pipeline_outputs
for sample in TL-22-JICYR8PP_N_DSQ1 TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1 TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
    if [ ! -e /Volumes/UTUCproject/DNA_Pipeline_outputs/${sample}/${sample}.dedup.numAligned.txt ]
    then
        cd ${sample}

        trim_galore --fastqc --paired -j 10 ${sample}_1.fastq ${sample}_2.fastq

        bwa mem -R '@RG\tID:1\tSM:'${sample}'\tPL:illumina\tLB:lib1\tPU:unit1' -t 15 /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa ${sample}_1_val_1.fq ${sample}_2_val_2.fq | samtools view -bS -F 0x900 - > ${sample}.aligned.bam
        samtools view -c -F 0x4 ${sample}.aligned.bam > ${sample}.numAligned.txt

        rm -R ${sample}_1_val_1.fq ${sample}_2_val_2.fq

        mkdir tmp
        samtools sort -@ 10 -m 5G -T tmp/aln.sorted -o ${sample}.aligned.bam ${sample}.aligned.bam
        samtools index ${sample}.aligned.bam

        picard MarkDuplicates --INPUT ${sample}.aligned.bam --OUTPUT ${sample}.dedup.aligned.bam --METRICS_FILE ${sample}.metrics.txt --ASSUME_SORTED True --REMOVE_DUPLICATES True
        samtools sort -@ 10 -m 3G -T tmp/aln.sorted -o ${sample}.dedup.aligned.bam ${sample}.dedup.aligned.bam
        samtools index ${sample}.dedup.aligned.bam

        samtools view -c -F 0x4 ${sample}.dedup.aligned.bam > ${sample}.dedup.numAligned.txt
        cd ..

    fi
done

for sample in TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1
do
    cd ${sample}
    gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I ${sample}.dedup.aligned.bam -I /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -normal TL-22-JICYR8PP_N_DSQ1 -O ${sample}.GATK.somatic.dedup.vcf
    grep -w '^#\|chr[1-9]\|chr[1-2][0-9]\|chr[X]' ${sample}.GATK.somatic.dedup.vcf > CLEAN_${sample}.GATK.somatic.dedup.vcf
    cd ..
done

for sample in TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
    cd ${sample}
    gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I ${sample}.dedup.aligned.bam -O ${sample}.GATK.somatic.dedup.vcf
    grep -w '^#\|chr[1-9]\|chr[1-2][0-9]\|chr[X]' ${sample}.GATK.somatic.dedup.vcf > CLEAN_${sample}.GATK.somatic.dedup.vcf
    cd ..
done

cd TL-22-JICYR8PP_N_DSQ1
gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.aligned.bam -O GERMLINE_TL-22-JICYR8PP_N.GATK.somatic.vcf
gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -O GERMLINE_TL-22-JICYR8PP_N.GATK.somatic.dedup.vcf
cd ..


#testing freebayes to run through COSMIC
for sample in TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1
do
    cd ${sample}
    freebayes -F 0.01 -C 2 --pooled-continuous --fasta-reference /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa ${sample}.aligned.bam /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.aligned.bam > freebayes_${sample}_variants.vcf
    freebayes -F 0.01 -C 2 --pooled-continuous --fasta-reference /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa ${sample}.dedup.aligned.bam /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam > freebayes_${sample}_variants.dedup.vcf
    cd ..
done

for sample in TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
    cd ${sample}
    freebayes -F 0.01 -C 2 --pooled-continuous --fasta-reference /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa ${sample}.aligned.bam > freebayes_${sample}_variants.vcf
    freebayes -F 0.01 -C 2 --pooled-continuous --fasta-reference /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa ${sample}.dedup.aligned.bam > freebayes_${sample}_variants.dedup.vcf
    cd ..
done



#run in main env

cd /Volumes/UTUCproject/DNA_Pipeline_outputs
sequenza-utils gc_wiggle -f /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -o sequenza_gc_wiggle.txt

cd /Volumes/UTUCproject/DNA_Pipeline_outputs
for sample in TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1
do
    cd ${sample}
    sequenza-utils bam2seqz -n /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -gc /Volumes/UTUCproject/DNA_Pipeline_outputs/sequenza_gc_wiggle.txt -t ${sample}.dedup.aligned.bam -F /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -o ${sample}.sequenza.txt
    cd ..
done

#run in feature env
configManta.py --normalBam=/Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.aligned.bam --tumorBam=/Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-86QRKCES_T_DSQ1/TL-22- 86QRKCES_T_DSQ1.aligned.bam --referenceFasta=/Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa