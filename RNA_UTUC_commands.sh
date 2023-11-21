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


for sample in $(ls -d */ | cut -f1 -d'/')
do
	if [ ! -d /temp_data/RNA_data/${sample}/RMATS/summary.txt ]
	then
		cd ${sample}	
		mkdir tmp
		mkdir RMATS
		#running it by itself first using --statoff but then remove that flag in order to compare within each patient 
		echo "/temp_data/RNA_data/${sample}/${sample}_1.fastq:/temp_data/RNA_data/${sample}/${sample}_3.fastq" >> s1.txt
		python /home/rin/anaconda3/envs/rmats_env/rMATS/rmats.py --s1 s1.txt --gtf /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --bi /home/rin/Desktop/Reference_files/hg19_files_goldenpath -t paired --readLength 50 --nthread 28 --od /temp_data/RNA_data/${sample}/RMATS --tmp /temp_data/RNA_data/${sample}/tmp --statoff
		cd ..
	fi
done















#SASHIMI

conda activate pipeline_macvs
cd /temp_data/RNA_data
ulimit -n 400000

echo "TL-22-JD5PK2C7_T_RSQ1	TL-22-JD5PK2C7_T_RSQ1/TL-22-JD5PK2C7_T_RSQ1_RNAAligned.sortedByCoord.out.bam" >> sashimi_groups_file.txt
echo "TL-23-FZU8JWAH_T_RSQ1	TL-23-FZU8JWAH_T_RSQ1/TL-23-FZU8JWAH_T_RSQ1_RNAAligned.sortedByCoord.out.bam" >> sashimi_groups_file.txt

/home/rin/ggsashimi.py -b sashimi_groups_file.txt -c chr13:32889607-32974509 -o SASHIMI_BRCA2whole -g /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --height 4 --width 20 -F pdf -C 2 -P palette.txt --shrink
/home/rin/ggsashimi.py -b sashimi_groups_file.txt -c chr13:32889607-32974509 -o SASHIMI_BRCA2whole -g /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --height 4 --width 20 -F jpeg -C 2 -P palette.txt --shrink

/home/rin/ggsashimi.py -b sashimi_groups_file.txt -c chr13:32919825-32935429 -o SASHIMI_BRCA2exons13_16 -g /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --height 4 --width 20 -F pdf -C 2 -P palette.txt --shrink
/home/rin/ggsashimi.py -b sashimi_groups_file.txt -c chr13:32919825-32935429 -o SASHIMI_BRCA2exons13_16 -g /home/rin/Desktop/Reference_files/hg19_files_goldenpath/hg19.refGene.gtf --height 4 --width 20 -F jpeg -C 2 -P palette.txt --shrink


