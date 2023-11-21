cd /media/rin/UTUCproject/RNA_Pipeline_outputs
ulimit -n 40000

conda activate pipeline_macvs

for sample in $(ls -d */ | cut -f1 -d'/')
do
if [ ! -e /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/RSEM_${sample}_RNAAligned.toTranscriptome.out.bam ]
then
	cd ${sample}
	STAR --genomeDir /media/rin/UTUCproject/RNA_Pipeline_outputs \
--runThreadN 26 \
--readFilesIn /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/${sample}_1.fastq /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/${sample}_3.fastq \
--outFileNamePrefix /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/${sample}_RNA \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard --quantMode GeneCounts
	samtools index /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/${sample}_RNAAligned.sortedByCoord.out.bam
	
	
	STAR --genomeDir /media/rin/UTUCproject/RNA_Pipeline_outputs \
--runThreadN 26 \
--readFilesIn /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/${sample}_1.fastq /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/${sample}_3.fastq \
--outFileNamePrefix /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/RSEM_${sample}_RNA \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard --quantMode TranscriptomeSAM
	
	rsem-calculate-expression --paired-end --alignments -p 26 /media/rin/UTUCproject/RNA_Pipeline_outputs/${sample}/RSEM_${sample}_RNAAligned.toTranscriptome.out.bam /media/rin/UTUCproject/RNA_Pipeline_outputs/hg19 ${sample}
	cd ..

fi
done






#DONT RUN CHECK FIRST IF YOU NEED IT
rsem-prepare-reference --gtf  /media/rin/UTUCproject/RNA_Pipeline_outputs/hg19.refGene.gtf \
    --star \
    --star-sjdboverhang 99 \
     -p 8 \
      /media/rin/UTUCproject/RNA_Pipeline_outputs/hg19.fa \
      /media/rin/UTUCproject/RNA_Pipeline_outputs/hg19
      

#STAR fusion
cd /home/rin/Desktop/UTUCproject/RNA_Pipeline_outputs
conda activate starfus_env
ulimit -n 400000

for sample in $(ls -d */ | cut -f1 -d'/')
do
	if [ ! -d /temp_data/${sample}/STARFUS_${sample} ]
	then
		cd ${sample}	
		STAR-Fusion --left_fq ${sample}_1.fastq  --right_fq ${sample}_3.fastq --genome_lib_dir /home/rin/Desktop/Reference_files/GRCh37_gencode_v19_CTAT_lib_Mar272019.plug-n-play/ctat_genome_lib_build_dir --CPU 30 --output_dir STARFUS_${sample}
		cd ..
	fi
done



