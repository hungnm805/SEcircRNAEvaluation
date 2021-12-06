#### Command line for execution of circRNA detection methods

### <I>. Alignment
## Some methods use some alignment setting from aligners like bwa or STAR.

##1. BWA
##version: 0.7.17-r1188
#Indexing genome

bwa index $FASTA_genome

##2. STAR
##version: 2.7.5c
# Indexing genome

STAR --runThreadN $num_of_threads --runMode genomeGenerate --genomeDir $STAR_index_dir --genomeFastaFiles --sjdbGTFfile $GTF_annotation_file --genomeSAsparseD 2 --genomeSAindexNbases 14

# Alignment (Single-end)

STAR --chimSegmentMin 20 --runThreadN $num_of_threads --genomeDir $STAR_index_dir --readFilesIn $FASTQ_seq_se --readFilesCommand zcat --chimScoreMin 1 --alignIntronMax 1000000 --outFilterMismatchNoverReadLmax 0.02 --alignTranscriptsPerReadNmax 100000 --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --outFilterMultimapNmax 2 --limitBAMsortRAM 22298667641 --outFileNamePrefix $STAR_output_path

# Alignment (Paired-End)

STAR --chimSegmentMin 20 --runThreadN $num_of_threads --genomeDir $STAR_index_dir --readFilesIn $FASTQ_seq_1 $FASTQ_seq_2 --readFilesCommand zcat --chimScoreMin 1 --alignIntronMax 1000000 --outFilterMismatchNoverReadLmax 0.02 --alignTranscriptsPerReadNmax 100000 --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --outFilterMultimapNmax 2 --limitBAMsortRAM 22298667641 --outFileNamePrefix $STAR_output_path

### <II>. circRNA detection

## 1. CIRI2
##version: 2.0.6
##URL: https://sourceforge.net/projects/ciri/

#Single-end

bwa mem -t $num_of_threads $FASTA_genome $FASTQ_seq_se \
| perl CIRI2.pl -I - -O $CR2_result_file -F $FASTA_genome -A $GTF_annotation_file -0 -T $num_of_threads

#Paired-End

bwa mem -t $num_of_threads $FASTA_genome $FASTQ_seq_1 $FASTQ_seq_2 \
| perl CIRI2.pl -I - -O $CR2_result_file -F $FASTA_genome -A $GTF_annotation_file -0 -T $num_of_threads

## 2. CIRCExplorer
##version: 1.1.10
##URL: https://github.com/YangLab/CIRCexplorer

star_parse.py $STAR_output_path/Chimeric.out.junction $STAR_output_path/fusion_junction.txt

CIRCexplorer.py -j $STAR_output_path/fusion_junction.txt -g -r $CE_ref_txt -o $CE_output_path

## 3. DCC
##version: 0.5.0
##URL: https://github.com/dieterich-lab/DCC

samtools index -b $STAR_output_path/Aligned.sortedByCoord.out.bam $STAR_output_path/Aligned.sortedByCoord.out.bam.bai

DCC -T $num_of_threads $STAR_output_path/Chimeric.out.junction -D -R $DCC_repeat_gtf -an $GTF_annotation_file -B $STAR_output_path/Aligned.sortedByCoord.out.bam -F -M -Nr 1 1 -fg -G -A -O $DCC_output_dir

## 4. CircRNA_finder
##version: 1.2
##URL: https://github.com/orzechoj/circRNA_finder

perl postProcessStarAlignment.pl --starDir $STAR_output_path --minLen 1 --outDir $CF_output_path

## 5. Find_circ2
##version: 1.99
##URL: https://github.com/rajewsky-lab/find_circ2

#Single-end

bwa mem -t $num_of_threads -A2 -B10 -k 15 -T 1 $FASTA_genome $FASTQ_seq_se \
| find_circ.py --genome $FASTA_genome -o $FC2_output_dir

#Paired-End

bwa mem -t$num_of_threads -A2 -B10 -k 15 -T 1 $FASTA_genome $FASTQ_seq_1 $FASTQ_seq_2 \
| find_circ.py --genome $FASTA_genome -o $FC2_output_dir

## 6. UROBORUS
##version: 2.0.0
##URL: https://github.com/WGLab/UROBORUS

#Single_end

tophat --bowtie1 -p $num_of_threads -o $UR_output_dir $bowtie1_index_dir $FASTQ_seq_se

samtools view $UR_output_dir/unmapped.bam > $UR_output_dir/unmapped.sam

perl UROBORUS.pl -p $num_of_threads -index $bowtie1_index_dir -gtf $GTF_annotation_file -fasta $chr_genome_dir $UR_output_dir/unmapped.sam $UR_output_dir/accepted_hits.bam

#Paired-end

tophat --bowtie1 -p $num_of_threads -o $UR_output_dir $bowtie1_index_dir $FASTQ_seq_1 $FASTQ_seq_2 \

samtools view $UR_output_dir/unmapped.bam > $UR_output_dir/unmapped.sam

perl UROBORUS.pl -p $num_of_threads -index $bowtie1_index_dir -gtf $GTF_annotation_file -fasta $chr_genome_dir $UR_output_dir/unmapped.sam $UR_output_dir/accepted_hits.bam

## 7. Clirc
##version: 0.1.0
##URL: https://github.com/Minzhe/Clirc

perl Clirc_library.pl -coord $circbase_circRNA_list -genome $FASTA_genome -library $Clirc_lib_dir

perl Clirc_search.pl -thread $num_of_threads -fastq $FASTQ_seq_se -library $Clirc_lib_dir -results $Clirc_result_file -cleanup

perl Clirc_filter.pl -out $Clirc_result_file -input $Clirc_result_file -library $Clirc_lib_dir


## 8. CircScan
##version: 0.1
##URL: https://github.com/sysu-software/circscan

bowtie2 -p $num_of_threads -t -k 4 --no-unal -D 200 -R 3 -N 0 -L 15 -i S,1,0.5 --score -min=C,-16,0 -q -U $FASTQ_seq_se -x $bowtie2_index_dir --un $bowtie2_alignment_path_minus.unmapped.fastq > $bowtie2_alignment_path_minus.alignments.bwt2 2> $bowtie2_alignment_path_minus.log

samtools view -bS $bowtie2_alignment_path_minus.alignments.bwt2 \
|samtools sort - \
|samtools index - > $bowtie2_alignment_path_minus.alignments.sorted.bam

bowtie2 -p $num_of_threads -t --no-unal -D 200 -R 3 -N 0 -L 15 -i S,1,0.5 --score -min=C,-16,0 -q -U $bowtie2_alignment_path_minus.unmapped.fastq -x <bowtie2 transcriptome index directory > --un $bowtie2_alignment_path_minus.transcriptome.unmapped.fastq $bowtie2_alignment_path_minus.transcriptome.alignments.bwt2

reformat.sh in=$bowtie2_alignment_path_minus.transcriptome.unmapped.fastq out=$bowtie2_alignment_path_minus.transcriptome.unmapped.fa overwrite=T

circScan $FASTA_genome $genome_index_fai $BED12_annotation $bowtie2_alignment_path_minus.alignments.sorted.bam $bowtie2_alignment_path_minus.transcriptome.unmapped.fa > $CS_result_file

