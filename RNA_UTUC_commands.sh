cd /temp_data/RNA_data
ulimit -n 400000
conda activate pipeline_macvs

for sample in $(ls -d */ | cut -f1 -d'/')
do
if [ ! -e /temp_data/RNA_data/${sample}/RSEM_${sample}_RNAAligned.toTranscriptome.out.bam ]
then
	cd ${sample}
	STAR --genomeDir /home/rin/Desktop/Reference_files/hg19_files_goldenpath \
--runThreadN 26 \
--readFilesIn /temp_data/RNA_data/${sample}/${sample}_1.fastq /temp_data/RNA_data/${sample}/${sample}_3.fastq \
--outFileNamePrefix /temp_data/RNA_data/${sample}/${sample}_RNA \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard --quantMode GeneCounts
	samtools index /temp_data/RNA_data/${sample}/${sample}_RNAAligned.sortedByCoord.out.bam
	
	
	STAR --genomeDir /home/rin/Desktop/Reference_files/hg19_files_goldenpath \
--runThreadN 26 \
--readFilesIn /temp_data/RNA_data/${sample}/${sample}_1.fastq /temp_data/RNA_data/${sample}/${sample}_3.fastq \
--outFileNamePrefix /temp_data/RNA_data/${sample}/RSEM_${sample}_RNA \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard --quantMode TranscriptomeSAM
	
	rsem-calculate-expression --paired-end --alignments -p 26 /temp_data/RNA_data/${sample}/RSEM_${sample}_RNAAligned.toTranscriptome.out.bam /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19 ${sample}
	cd ..

fi
done

#TOPHAT
cd /temp_data/RNA_data
conda activate tophat_env
ulimit -n 400000

for sample in $(ls -d */ | cut -f1 -d'/')
do
	if [ ! -e /temp_data/RNA_data/${sample}/tophat_${sample}/fusions.out ]
	then
		cd ${sample}	
		tophat -o tophat_${sample} -p 30 --fusion-search --keep-fasta-order --no-coverage-search -r 0 --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 /home/rin/Desktop/Reference_files/tophat-directory/hg19 ${sample}_1.fastq ${sample}_3.fastq
		tophat-fusion-post -p 30 -o tophat_${sample}/ /home/rin/Desktop/Reference_files/tophat-directory/hg19
		cd ..
	fi
done


#STARFUSION
cd /temp_data/RNA_data
conda activate starfus_env
ulimit -n 400000

for sample in $(ls -d */ | cut -f1 -d'/')
do
	if [ ! -d /temp_data/RNA_data/${sample}/STARFUS_${sample} ]
	then
		cd ${sample}	
		STAR-Fusion --left_fq ${sample}_1.fastq  --right_fq ${sample}_3.fastq --genome_lib_dir /home/rin/Desktop/Reference_files/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir --CPU 26 --output_dir STARFUS_${sample}
		cd ..
	fi
done


#RMATS
cd /temp_data/RNA_data
conda activate rmats_env
ulimit -n 400000


for sample in TL-22-354IF4A9_T_RSQ1 TL-22-86QRKCES_T_RSQ1 TL-22-BK3DDUYG_T_RSQ1 TL-22-DHERTUS6_T_RSQ1 TL-22-JICYR8PP_T_RSQ1 TL-22-KXMGA3WE_T_RSQ1 TL-22-P4ACNIE3_T_RSQ1 TL-22-PKA8ZUD2_T_RSQ1 TL-22-RVHYDX4F_T_RSQ1 TL-22-SA4HH23W_T_RSQ1 TL-22-TMY85WNT_T_RSQ1
do
	if [ ! -d /temp_data/RNA_data/${sample}/RMATS/summary.txt ]
	then
		cd ${sample}	
		mkdir RMATS
		mkdir RMATS/tmp
		#running it by itself first using --statoff but then remove that flag in order to compare within each patient 
		echo "/temp_data/RNA_data/${sample}/${sample}_1.fastq:/temp_data/RNA_data/${sample}/${sample}_3.fastq" > s1.txt
		seqLength=$(head -n 2 /temp_data/RNA_data/${sample}/${sample}_1.fastq | tail -n 1 | wc -c)
		python /home/rin/anaconda3/envs/rmats_env/rMATS/rmats.py --s1 s1.txt --gtf /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --bi /home/rin/Desktop/Reference_files/hg19_files_goldenpath -t paired --readLength 150 --variable-read-length --nthread 28 --od /temp_data/RNA_data/${sample}/RMATS --tmp /temp_data/RNA_data/${sample}/RMATS/tmp --statoff
		cd ..
	fi
done


#RMATS comparing each sample to each other
cd /temp_data/RNA_data
conda activate rmats_env
ulimit -n 400000
for sample_1 in TL-22-354IF4A9_T_RSQ1 TL-22-86QRKCES_T_RSQ1 TL-22-BK3DDUYG_T_RSQ1 TL-22-DHERTUS6_T_RSQ1 TL-22-JICYR8PP_T_RSQ1 TL-22-KXMGA3WE_T_RSQ1 TL-22-P4ACNIE3_T_RSQ1 TL-22-PKA8ZUD2_T_RSQ1 TL-22-RVHYDX4F_T_RSQ1 TL-22-SA4HH23W_T_RSQ1 TL-22-TMY85WNT_T_RSQ1
do
	for sample_2 in TL-22-354IF4A9_T_RSQ1 TL-22-86QRKCES_T_RSQ1 TL-22-BK3DDUYG_T_RSQ1 TL-22-DHERTUS6_T_RSQ1 TL-22-JICYR8PP_T_RSQ1 TL-22-KXMGA3WE_T_RSQ1 TL-22-P4ACNIE3_T_RSQ1 TL-22-PKA8ZUD2_T_RSQ1 TL-22-RVHYDX4F_T_RSQ1 TL-22-SA4HH23W_T_RSQ1 TL-22-TMY85WNT_T_RSQ1
	do
	if [ ${sample_1} != ${sample_2} ]
	then
		if [ ! -d /temp_data/RNA_data/comparison_of_${sample_1}_and_${sample_2} ]
		then
			echo "sample1: ${sample_1}"
			echo "sample2: ${sample_2}"
			echo "/temp_data/RNA_data/${sample_1}/${sample_1}_1.fastq:/temp_data/RNA_data/${sample_1}/${sample_1}_3.fastq" > ${sample_1}_s1.txt
			echo "/temp_data/RNA_data/${sample_2}/${sample_2}_1.fastq:/temp_data/RNA_data/${sample_2}/${sample_2}_3.fastq" > ${sample_2}_s2.txt
			mkdir comparison_of_${sample_1}_and_${sample_2}
			mkdir comparison_of_${sample_1}_and_${sample_2}/tmp
			seqLength=$(head -n 2 /temp_data/RNA_data/${sample_1}/${sample_1}_1.fastq | tail -n 1 | wc -c)
			python /home/rin/anaconda3/envs/rmats_env/rMATS/rmats.py --s1 ${sample_1}_s1.txt --s2 ${sample_2}_s2.txt --gtf /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --bi /home/rin/Desktop/Reference_files/hg19_files_goldenpath -t paired --readLength 150 --variable-read-length --nthread 28 --od /temp_data/RNA_data/comparison_of_${sample_1}_and_${sample_2} --tmp /temp_data/RNA_data/comparison_of_${sample_1}_and_${sample_2}/tmp 
		fi
	
	fi
	done
done



#SASHIMI
conda activate pipeline_macvs
cd /temp_data/RNA_data
ulimit -n 400000

for sample in TL-22-TMY85WNT_T_RSQ1 TL-22-KXMGA3WE_T_RSQ1 TL-22-BK3DDUYG_T_RSQ1
do
	echo "${sample}	${sample}/${sample}_RNAAligned.sortedByCoord.out.bam" >> sashimi_groups_file_FGFR3TACC3.txt
done

/home/rin/ggsashimi.py -b sashimi_groups_file_FGFR3TACC3.txt -c chr4:1741429-1808661 -o SASHIMI_FGFR3_TACC3 -g /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --height 4 --width 20 -F pdf --shrink
/home/rin/ggsashimi.py -b sashimi_groups_file_FGFR3TACC3.txt -c chr4:1741429-1808661 -o SASHIMI_FGFR3_TACC3 -g /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --height 4 --width 20 -F jpeg --shrink


conda activate pipeline_macvs
cd /temp_data/RNA_data
ulimit -n 400000

for sample in TL-22-TMY85WNT_T_RSQ1 TL-22-KXMGA3WE_T_RSQ1 TL-22-BK3DDUYG_T_RSQ1
do
	python /home/rin/rmats2sashimiplot/src/rmats2sashimiplot/rmats2sashimiplot.py --b1 ${sample}/${sample}_RNAAligned.sortedByCoord.out.bam --b2 ${sample}/${sample}_RNAAligned.sortedByCoord.out.bam -c chr16:+:9000:25000:/home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --l1 SampleOne --l2 SampleTwo --exon_s 1 --intron_s 5 -o test_coordinate_output
done








