#!/bin/bash

# this script is used to create a bunch of files estimating pcr / sequencing error rates in the real data. 
# it also creates dataframes of the vcf files, which I find are a bit easier to work with than the vcf files. 
# you might need to install the packages jupyter nbconvert (and samtools and bcftools and bowtie2)

shopt -s extglob
#FILES="Sample_Ldn_110121_rep1 Sample_Ldn_110121_rep2 Sample_Ldn_120121_rep1 Sample_Ldn_120121_rep2 Sample_Ldn_130121_rep1 Sample_Ldn_130121_rep2 Sample_Ldn_140121_rep1 Sample_Ldn_140121_rep2"
#FILES="Sample_Ldn_110121_rep1 Sample_Ldn_110121_rep2 Sample_Ldn_120121_rep1 Sample_Ldn_120121_rep2 Sample_Ldn_130121_rep1 Sample_Ldn_130121_rep2 Sample_Ldn_140121_rep1 Sample_Ldn_140121_rep2 Sample_1_211283_A1 Sample_1_211283_A2 Sample_1_211283_A3 Sample_1_211283_A4 Sample_1_211283_A5 Sample_1_211283_A6 Sample_1_211283_A7 Sample_1_211283_A8 Sample_1_211283_A9 Sample_1_211283_A10 Sample_1_211283_A11 Sample_1_211283_A12 Sample_1_211283_B1 Sample_1_211283_B2 Sample_1_211283_B3 Sample_1_211283_B4 Sample_1_211283_B5 Sample_1_211283_B6 Sample_1_211283_B7 Sample_1_211283_B8 Sample_1_211283_B9 Sample_1_211283_B10 Sample_1_211283_B11 Sample_1_211283_B12 Sample_1_211283_C1 Sample_1_211283_C2 Sample_1_211283_C3 Sample_1_211283_C4 Sample_1_211283_C5 Sample_1_211283_C6 Sample_1_211283_C7 Sample_1_211283_C8 Sample_1_211283_C9 Sample_1_211283_C10 Sample_1_211283_C11 Sample_1_211283_C12 Sample_1_211283_D1 Sample_1_211283_D2 Sample_1_211283_D3 Sample_1_211283_D4 Sample_1_211283_D5 Sample_1_211283_D6 Sample_1_211283_D7 Sample_1_211283_D8 Sample_1_211283_D9 Sample_1_211283_D10 Sample_1_211283_D11 Sample_1_211283_D12 Sample_1_211283_E1 Sample_1_211283_E2 Sample_1_211283_E3 Sample_1_211283_E4 Sample_1_211283_E5 Sample_1_211283_E6 Sample_1_211283_E7 Sample_1_211283_E8 Sample_1_211283_E9 Sample_1_211283_E10 Sample_1_211283_E11 Sample_1_211283_E12 Sample_1_211283_F1 Sample_1_211283_F2 Sample_1_211283_F3 Sample_1_211283_F4 Sample_1_211283_F5 Sample_1_211283_F6 Sample_1_211283_F7 Sample_1_211283_F8 Sample_1_211283_F9 Sample_1_211283_F10 Sample_1_211283_F11 Sample_1_211283_F12 Sample_1_211283_G1 Sample_1_211283_G2 Sample_1_211283_G3 Sample_1_211283_G4 Sample_1_211283_G5 Sample_1_211283_G6 Sample_1_211283_G7 Sample_1_211283_G8 Sample_1_211283_G9 Sample_1_211283_G10 Sample_1_211283_G11 Sample_1_211283_G12 Sample_1_211283_H1 Sample_1_211283_H2 Sample_1_211283_H3 Sample_1_211283_H4 Sample_1_211283_H5 Sample_1_211283_H6 Sample_1_211283_H7 Sample_1_211283_H8 Sample_1_211283_H9 Sample_1_211283_H10 Sample_1_Control_Negative_H11 Sample_1_Control_Positive_H12 Sample_1_211262_A1 Sample_1_211262_A2 Sample_1_211262_A3 Sample_1_211262_A4 Sample_1_211262_A5 Sample_1_211262_A6 Sample_1_211262_A7 Sample_1_211262_A8 Sample_1_211262_A9 Sample_1_211262_A10 Sample_1_211262_B1 Sample_1_211262_B2 Sample_1_211262_B3 Sample_1_211262_B4 Sample_1_211262_B5 Sample_1_211262_B6 Sample_1_211262_B7 Sample_1_211262_B8 Sample_1_211262_B9 Sample_1_211262_B10 Sample_1_211262_B11 Sample_1_211262_B12 Sample_1_211262_C1 Sample_1_211262_C2 Sample_1_211262_C3 Sample_1_211262_C4 Sample_1_211262_C5 Sample_1_211262_C6 Sample_1_211262_C7 Sample_1_211262_C8 Sample_1_211262_C9 Sample_1_211262_C10 Sample_1_211262_C11 Sample_1_211262_C12 Sample_1_211262_D1 Sample_1_211262_D2 Sample_1_211262_D3 Sample_1_211262_D4 Sample_1_211262_D5 Sample_1_211262_D6 Sample_1_211262_D7 Sample_1_211262_D8 Sample_1_211262_D9 Sample_1_211262_D10 Sample_1_211262_D11 Sample_1_211262_D12 Sample_1_211262_E1 Sample_1_211262_E2 Sample_1_211262_E3 Sample_1_211262_E4 Sample_1_211262_E5 Sample_1_211262_E6 Sample_1_211262_E7 Sample_1_211262_E8 Sample_1_211262_E9 Sample_1_211262_E10 Sample_1_211262_E11 Sample_1_211262_E12 Sample_1_211262_F1 Sample_1_211262_F2 Sample_1_211262_F3 Sample_1_211262_F4 Sample_1_211262_F5 Sample_1_211262_F6 Sample_1_211262_F7 Sample_1_211262_F8 Sample_1_211262_F9 Sample_1_211262_F10 Sample_1_211262_F11 Sample_1_211262_F12 Sample_1_211262_G1 Sample_1_211262_G2 Sample_1_211262_G3 Sample_1_211262_G4 Sample_1_211262_G5 Sample_1_211262_G6 Sample_1_211262_G7 Sample_1_211262_G8 Sample_1_211262_G9 Sample_1_211262_G10 Sample_1_211262_G11 Sample_1_211262_G12 Sample_1_211262_H1 Sample_1_211262_H2 Sample_1_211262_H3 Sample_1_211262_H4 Sample_1_211262_H5 Sample_1_211262_H6 Sample_1_211262_H7 Sample_1_211262_H8 Sample_1_211262_H9 Sample_1_211262_H10 Sample_1_211262_H11 Sample_1_211262_H12 Sample_1_Control_Negative_A11 Sample_1_Control_Positive_A12"
FILES="A100_rep2_highconc Sample_1_211283_F9 Sample_1_211283_F1 Sample_1_211262_C8 Sample_1_211283_A12 Sample_1_211262_C9 Sample_1_211262_E12 A0_rep2_lowconc Sample_1_211262_G3 Sample_1_211283_A1 Sample_1_211262_A3 Sample_1_211262_E11 Sample_1_211283_H1 Sample_1_211283_C8 Sample_1_211262_E10 Sample_1_211262_A1 Sample_1_211283_H9 Sample_Ldn_110121_rep1 Sample_1_211262_D2 Sample_1_211283_F8 Sample_1_211283_G1 Sample_1_211262_E4 Sample_1_211283_F5 Sample_1_211262_H11 Sample_1_211283_D1 Sample_Ldn_120121_rep2 Sample_1_211283_C1 Sample_1_211262_F9 Sample_1_211283_E2 Sample_1_211262_D5 Sample_1_211283_D7 Sample_1_211283_B4 Sample_1_211283_F4 Sample_1_211262_H4 Sample_1_211262_D4 Sample_Ldn_110121_rep2 Sample_1_211283_E5 Sample_1_211262_G10 Sample_1_211283_C7 Sample_1_211262_H8 Sample_1_211283_A8 Sample_1_211283_C6 Sample_1_211283_A11 Sample_Ldn_120121_rep1 Sample_1_211262_D8 Sample_1_211283_G3 Sample_1_211262_F1 Sample_1_211283_G8 Sample_1_Control_Positive_A12 Sample_1_Control_Negative_A11 Sample_1_211283_H8 Sample_1_211283_G10 A100_rep2_lowconc Sample_1_211262_A2 Sample_1_211262_E5 Sample_1_Control_Positive_H12 Sample_1_211262_G9 Sample_1_211283_C2 Sample_1_211262_F2 Sample_1_211283_A9 Sample_1_211262_F8 Sample_1_211262_A4 Sample_1_211262_E3 Sample_1_211262_C12 Sample_1_211262_C5 Sample_1_211283_C12 Sample_1_211262_B8 Sample_1_211262_B11 Sample_1_211262_E7 Sample_1_211283_C9 Sample_1_211283_B7 Sample_1_211262_A8 Sample_1_211283_G6 Sample_1_211262_D7 Sample_1_211283_G12 Sample_1_211283_E8 Sample_1_211262_H5 Sample_1_211283_E10 Sample_1_211283_C5 Sample_1_211262_G6 Sample_1_211283_A5 Sample_1_211283_H3 Sample_1_211262_F10 Sample_1_211262_E9 Sample_1_211283_G2 Sample_1_211283_B8 Sample_1_211262_C3 Sample_1_211262_G12 A1_rep1_highconc Sample_1_211283_B12 Sample_1_211262_G2 Sample_1_211262_E1 Sample_1_211283_D5 A0_rep1_highconc Sample_1_211283_F12 Sample_1_211283_A4 Sample_1_211262_C6 Sample_Ldn_130121_rep1 Sample_1_211283_D9 Sample_1_211262_H12 Sample_1_211262_B9 Sample_1_211262_B5 Sample_1_211262_D11 Sample_1_211262_B10 Sample_1_211262_D1 Sample_1_211283_C10 Sample_1_211262_H3 Sample_1_211283_F6 Sample_1_211262_A5 Sample_1_211283_E12 Sample_1_211262_G8 Sample_1_211283_E9 Sample_1_211283_B10 Sample_1_211283_F11 Sample_1_211262_G11 Sample_1_211283_H7 Sample_1_211262_H2 Sample_1_211283_B1 A0_rep2_highconc Sample_1_211283_A6 Sample_1_211283_G4 Sample_1_211283_F3 Sample_1_211283_B3 Sample_1_211262_F7 Sample_1_211262_C10 Sample_1_211262_H1 Sample_1_211262_G1 Sample_1_211283_E3 Sample_1_211283_E6 Sample_1_211262_E8 Sample_1_211283_E11 Sample_1_211262_F12 Sample_1_211262_A9 Sample_1_211262_C4 Sample_1_211262_D6 Sample_1_211262_A7 Sample_1_211262_C7 Sample_1_Control_Negative_H11 A1_rep2_highconc Sample_1_211262_E2 Sample_1_211283_G9 Sample_1_211262_D9 Sample_1_211283_A2 Sample_1_211283_D11 Sample_1_211262_D10 Sample_Ldn_140121_rep2 Sample_1_211283_C4 Sample_1_211262_G7 Sample_1_211283_D6 Sample_1_211283_C3 Sample_1_211283_H5 Sample_1_211283_H2 Sample_1_211262_B6 Sample_1_211283_F7 Sample_1_211283_D2 Sample_1_211262_G5 Sample_1_211283_F2 Sample_1_211262_F11 Sample_1_211262_F3 Sample_1_211283_A3 Sample_1_211283_H6 Sample_1_211283_E7 Sample_1_211262_F4 Sample_1_211262_E6 Sample_1_211283_B5 Sample_1_211283_D12 Sample_1_211283_B11 Sample_1_211262_H9 Sample_1_211262_C1 A0_rep1_lowconc Sample_1_211262_C11 Sample_Ldn_100121_rep2 Sample_Ldn_130121_rep2 Sample_1_211283_H4 Sample_1_211283_G11 Sample_1_211283_D4 Sample_1_211262_G4 Sample_1_211283_H10 Sample_1_211262_H6 Sample_Ldn_100121_rep1 A100_rep1_highconc Sample_1_211283_A10 Sample_1_211262_H7 Sample_1_211262_A6 Sample_1_211283_D8 Sample_1_211283_B6 Sample_1_211283_E1 Sample_1_211283_D10 Sample_1_211283_A7 Sample_1_211262_B1 Sample_1_211283_E4 Sample_1_211262_B2 Sample_1_211262_D3 Sample_1_211262_B12 Sample_1_211262_D12 Sample_1_211283_G7 Sample_1_211283_B2 Sample_1_211283_D3 Sample_1_211262_A10 Sample_1_211262_B3 Sample_Ldn_140121_rep1 Sample_1_211262_B4 A100_rep1_lowconc Sample_1_211283_B9 Sample_1_211283_F10 Sample_1_211262_B7 Sample_1_211262_C2 Sample_1_211262_H10 Sample_1_211283_G5 Sample_1_211262_F6 Sample_1_211283_C11 Sample_1_211262_F5"
#FILES="Sample_Ldn_100121_rep1 Sample_Ldn_100121_rep2"
#FILES="Sample_1_211283_A1 Sample_1_211283_A2"
for file in $FILES; do
    path=../data/ww_dump/plate_validation/$file
    reads=$(ls $path)
    set -- $reads
    
    # this code block might be necessary for some files - to move things from the backup folders to one level up
    #gzip -d $path/$1
    #gzip -d $path/$2
    #mv $path/backup/$1 $path/
    #mv $path/backup/$2 $path/
    #rmdir $path/backup
    
    bowtie2 -x ../data/indices/wuhan_indices/MN908947.3 -1 $path/$1 -2 $path/$2 | samtools view -bS - > $path/$file.bam
    samtools sort $path/$file.bam -o $path/$file.sorted.bam
    bcftools mpileup -f ../data/MN908947.3.fasta $path/$file.sorted.bam --annotate INFO/AD --max-depth 20000 -o $path/$file.vcf 
    
    #| bcftools view -Ov - > $path/$file_c15.bcf
    
    export FILENAME=$file
    jupyter nbconvert --to notebook --execute --inplace ../notes/frequency_spectrum_cleaned.ipynb
    unset FILENAME
done