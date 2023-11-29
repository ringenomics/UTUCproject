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

#SVVIZ
cd /temp_data/DNA_data
conda activate svviz_env
for sample in $(ls -d */ | cut -f1 -d'/')
do
    if [ ! -e /temp_data/DNA_data/${sample}/hi ]
    then
        cd ${sample}
        	svviz --type inv --processes 30 --annotations /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.bed -b ${sample}.dedup.aligned.bam /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa chr4 1740000 1810000
        cd ..
    fi
done

#FREEBAYEs on everyone
conda activate gatk_env
cd /intelpool/UTUCproject/DNA_data

for sample in TL-22-P4ACNIE3_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-86QRKCES_T_DSQ1 TL-22-354IF4A9_T_DSQ1 TL-22-RVHYDX4F_T_DSQ1 TL-22-SKGMA3U8_T_DSQ1 TL-22-KXMGA3WE_T_DSQ1 TL-22-BK3DDUYG_T_DSQ1 TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
    if [ ! -e /intelpool/UTUCproject/DNA_data/${sample}/filtered_freebayes_${sample}.bed ]
    then
        cd ${sample}
        	freebayes -f /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa ${sample}.dedup.aligned.bam > freebayes_${sample}.vcf
        	bcftools view -i 'QUAL>30 && FORMAT/DP>100' freebayes_${sample}.vcf > filtered_freebayes_${sample}.vcf
        	bcftools query -f '%CHROM\t%POS0\t%POS\t%ALT\n' filtered_freebayes_${sample}.vcf > filtered_freebayes_${sample}.bed
        cd ..
        
        cat ${sample}/filtered_freebayes_${sample}.bed >> freebayes_presorted_ALL_SAMPLES.bed
    fi
done

#now compare each sample within each patient
#pt1
for sample in TL-22-P4ACNIE3_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-86QRKCES_T_DSQ1
do
	cat ${sample}/filtered_freebayes_${sample}.bed >> pt1_freebayes_presorted_ALL_SAMPLES.bed
done


#pt2
for sample in TL-22-354IF4A9_T_DSQ1 TL-22-RVHYDX4F_T_DSQ1
do
	cat ${sample}/filtered_freebayes_${sample}.bed >> pt2_freebayes_presorted_ALL_SAMPLES.bed
done


#pt3
for sample in TL-22-SKGMA3U8_T_DSQ1 TL-22-KXMGA3WE_T_DSQ1 TL-22-BK3DDUYG_T_DSQ1
do
	cat ${sample}/filtered_freebayes_${sample}.bed >> pt3_freebayes_presorted_ALL_SAMPLES.bed
done


#pt4
for sample in TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
	cat ${sample}/filtered_freebayes_${sample}.bed >> pt4_freebayes_presorted_ALL_SAMPLES.bed
done

#running gatk now only on one patient with normal samples
conda activate gatk_env
cd /temp_data/DNA_data
for tumor_name in TL-22-P4ACNIE3_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-86QRKCES_T_DSQ1
do
	if [ ! -e /temp_data/DNA_data/${tumor_name}/final_${tumor_name}.GATK.somatic_filtered.dedup.bed ]
    	then
    		cd ${tumor_name}
		gatk Mutect2 -R /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa -I ${tumor_name}.dedup.aligned.bam -I /temp_data/DNA_data/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -normal TL-22-JICYR8PP_N_DSQ1 -O /temp_data/DNA_data/${tumor_name}/${tumor_name}.GATK.somatic_unfiltered.dedup.vcf --native-pair-hmm-threads 30
		gatk FilterMutectCalls -R /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa -V /temp_data/DNA_data/${tumor_name}/${tumor_name}.GATK.somatic_unfiltered.dedup.vcf -O /temp_data/DNA_data/${tumor_name}/${tumor_name}.GATK.somatic_filtered.dedup.vcf
		bcftools view -f 'PASS,.' ${tumor_name}.GATK.somatic_filtered.dedup.vcf > final_${tumor_name}.GATK.somatic_filtered.dedup.vcf
		
		bcftools query -f '%CHROM\t%POS0\t%POS\n' final_${tumor_name}.GATK.somatic_filtered.dedup.vcf > final_${tumor_name}.GATK.somatic_filtered.dedup.bed
		awk -v OFS="\t" '!/##/ {$10=$11=""}1' final_${tumor_name}.GATK.somatic_filtered.dedup.vcf |sed 's/^\s\+//g' > shortened_${tumor_name}.GATK.somatic_filtered.dedup.vcf
		
		cd ..
		
		cat ${tumor_name}/final_${tumor_name}.GATK.somatic_filtered.dedup.bed >> presorted_ALL_SAMPLES.bed

	fi
done

sort -k 1,1 -k2,2n presorted_ALL_SAMPLES.bed > sorted_ALL_SAMPLES.bed
bedtools merge -i sorted_ALL_SAMPLES.bed > mergedandsorted_ALL_SAMPLES.bed

#running pindel
conda activate gatk_env
cd /temp_data/DNA_data
for tumor_name in TL-22-P4ACNIE3_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-86QRKCES_T_DSQ1
do
	if [ ! -e /temp_data/DNA_data/${tumor_name}/FILTERED_pindel_output_${tumor_name}.vcf ]
    	then
		cd ${tumor_name}
		#making necessary pindel files
		echo "/temp_data/DNA_data/${tumor_name}/${tumor_name}.dedup.aligned.bam	250	${tumor_name}" >> ${tumor_name}_pindelinput.txt
		echo "/temp_data/DNA_data/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam	250	TL-22-JICYR8PP_N_DSQ1" >> ${tumor_name}_pindelinput.txt
		pindel -f /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa -i ${tumor_name}_pindelinput.txt -o pindel_output_${tumor_name} -T 30 -w 1000 -c ALL
		pindel2vcf -P pindel_output_${tumor_name} -r /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa -R UCSCgoldenpathhg19 -d 20180821  -w 1000
		#filtering vcf files
		gatk IndexFeatureFile -I pindel_output_${tumor_name}.vcf
		grep "\(0\/1\|1\/1\|1\/0\|#\)" pindel_output_${tumor_name}.vcf > FILTERED_pindel_output_${tumor_name}.vcf
		cd ..
	fi
done





#FACETS CNV CALLS
conda activate pipeline_macvs
cd /temp_data/DNA_data
for tumor_name in TL-22-P4ACNIE3_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-86QRKCES_T_DSQ1
do
	if [ ! -e /temp_data/DNA_data/${tumor_name}/FACETS_CNV/${tumor_name} ]
    	then
		cd ${tumor_name}
		mkdir FACETS_CNV
		cnv_facets.R -t ${tumor_name}.aligned.bam -n /temp_data/DNA_data/TL-22-JICYR8PP_N_DSQ1/TL-22-JICYR8PP_N_DSQ1.dedup.aligned.bam -vcf /home/rin/Desktop/Reference_files/hg19_files_goldenpath/correctly_labelled_00-common_all.vcf.gz -o FACETS_CNV/${tumor_name} --snp-nprocs 30
		cd ..
	fi
done


#SOMALIER
conda activate gatk_env
cd /intelpool/UTUCproject/DNA_data

for sample in TL-22-P4ACNIE3_T_DSQ1 TL-22-JICYR8PP_T_DSQ1 TL-22-DHERTUS6_T_DSQ1 TL-22-86QRKCES_T_DSQ1 TL-22-354IF4A9_T_DSQ1 TL-22-RVHYDX4F_T_DSQ1 TL-22-SKGMA3U8_T_DSQ1 TL-22-KXMGA3WE_T_DSQ1 TL-22-BK3DDUYG_T_DSQ1 TL-22-PKA8ZUD2_T_DSQ1 TL-22-SA4HH23W_T_DSQ1
do
    if [ ! -e /intelpool/UTUCproject/DNA_data/somalier_cohort/hi ]
    then
        cd ${sample}
        	somalier extract -d /intelpool/UTUCproject/DNA_data/somalier_cohort/ --sites /home/rin/Desktop/Reference_files/hg19_files_goldenpath/somalier_sites.hg19.vcf -f /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.fa ${sample}.dedup.aligned.bam
        cd ..
    fi
done

somalier relate --infer /intelpool/UTUCproject/DNA_data/somalier_cohort/*.somalier


