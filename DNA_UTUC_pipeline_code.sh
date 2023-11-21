#Run in pipeline env

cd /home/rin/Desktop/UTUCproject
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


#running gatk
conda activate gatk_env
cd /home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs
for sample in TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1
do
        cd ${sample}
        gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I ${sample}.dedup.aligned.bam -I /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -normal TL-22-JICYR8PP_N_DSQ1 -O ${sample}.GATK.somatic.dedup.vcf --native-pair-hmm-threads 30
        gatk FilterMutectCalls -R /home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs/hg19.fa -V ${sample}.GATK.somatic.dedup.vcf -O ${sample}.filtered_GATK.somatic.dedup.vcf
        mkdir /home/rin/Desktop/UTUCproject/Data_Analysis/${sample}
        bcftools filter -O z -o /home/rin/Desktop/UTUCproject/Data_Analysis/${sample}/final_${sample}.GATK.somatic_filtered.dedup.vcf -i '%FILTER="PASS"' ${sample}.filtered_GATK.somatic.dedup.vcf
        cd ..
done

cd /home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs
for sample in TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
     cd ${sample}
     gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I ${sample}.dedup.aligned.bam -O ${sample}.GATK.somatic.dedup.vcf
     gatk FilterMutectCalls -R /home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs/hg19.fa -V ${sample}.GATK.somatic.dedup.vcf -O ${sample}.filtered_GATK.somatic.dedup.vcf
     mkdir /home/rin/Desktop/UTUCproject/Data_Analysis/${sample}
     bcftools filter -O z -o /home/rin/Desktop/UTUCproject/Data_Analysis/${sample}/final_${sample}.GATK.somatic_filtered.dedup.vcf -i '%FILTER="PASS"' ${sample}.filtered_GATK.somatic.dedup.vcf
     cd ..
done


if [ ! -e /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/GERMLINE_TL-22-JICYR8PP_N.GATK.somatic.dedup.vcf ]
then
    cd TL-22-JICYR8PP_N_DSQ1
    gatk Mutect2 -R /Volumes/UTUCproject/DNA_Pipeline_outputs/hg19.fa -I /Volumes/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -O GERMLINE_TL-22-JICYR8PP_N.GATK.somatic.dedup.vcf
    cd ..
fi


#manta run
conda activate strelka
cd /home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs
for sample in TL-22-86QRKCES_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-P4ACNIE3_T_DSQ1
do
        cd ${sample}
        /home/rin/anaconda3/envs/strelka/share/manta-1.6.0-2/bin/configManta.py --normalBam=/home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam --tumorBam=${sample}.dedup.aligned.bam --referenceFasta=/home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs/hg19.fa
        /home/rin/Desktop/UTUCproject/DNA_Pipeline_outputs/${sample}/MantaWorkflow/runWorkflow.py -m local -j 20
        cd ..
done
 

