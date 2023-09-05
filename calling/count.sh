#!/bin/bash

bedtools makewindows -g <(head -n 23 reference/human/GRCh38.size) -w 350 > sliding_window.bed

bedtools intersect -v -a sliding_window.bed -b reference/human/hg38-blacklist.v2.bed > sliding_window.tmp
mv sliding_window.tmp sliding_window.bed


mkdir bed

samtools view ../STARR/AA_D_1/AA_D_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AA_D_1.bed
samtools view ../STARR/AA_D_2/AA_D_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AA_D_2.bed
samtools view ../STARR/AA_R_1/AA_R_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AA_R_1.bed
samtools view ../STARR/AA_R_2/AA_R_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AA_R_2.bed

samtools view ../STARR/AK_D_1/AK_D_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AK_D_1.bed
samtools view ../STARR/AK_D_2/AK_D_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AK_D_2.bed
samtools view ../STARR/AK_R_1/AK_R_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AK_R_1.bed
samtools view ../STARR/AK_R_2/AK_R_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/AK_R_2.bed

samtools view ../STARR/KA_D_1/KA_D_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KA_D_1.bed
samtools view ../STARR/KA_D_2/KA_D_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KA_D_2.bed
samtools view ../STARR/KA_R_1/KA_R_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KA_R_1.bed
samtools view ../STARR/KA_R_2/KA_R_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KA_R_2.bed

samtools view ../STARR/KK_D_1/KK_D_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KK_D_1.bed
samtools view ../STARR/KK_D_2/KK_D_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KK_D_2.bed
samtools view ../STARR/KK_R_1/KK_R_1.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KK_R_1.bed
samtools view ../STARR/KK_R_2/KK_R_2.filtered.bam | awk '{OFS="\t"; if ($9>0) print $3,$4-1,$4-1+$9;}' > bed/KK_R_2.bed


mkdir count

bedtools coverage -counts -a sliding_window.bed -b bed/AA_D_1.bed > count/AA_D_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/AA_D_2.bed > count/AA_D_2.bed
bedtools coverage -counts -a sliding_window.bed -b bed/AA_R_1.bed > count/AA_R_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/AA_R_2.bed > count/AA_R_2.bed

bedtools coverage -counts -a sliding_window.bed -b bed/AK_D_1.bed > count/AK_D_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/AK_D_2.bed > count/AK_D_2.bed
bedtools coverage -counts -a sliding_window.bed -b bed/AK_R_1.bed > count/AK_R_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/AK_R_2.bed > count/AK_R_2.bed

bedtools coverage -counts -a sliding_window.bed -b bed/KA_D_1.bed > count/KA_D_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/KA_D_2.bed > count/KA_D_2.bed
bedtools coverage -counts -a sliding_window.bed -b bed/KA_R_1.bed > count/KA_R_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/KA_R_2.bed > count/KA_R_2.bed

bedtools coverage -counts -a sliding_window.bed -b bed/KK_D_1.bed > count/KK_D_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/KK_D_2.bed > count/KK_D_2.bed
bedtools coverage -counts -a sliding_window.bed -b bed/KK_R_1.bed > count/KK_R_1.bed
bedtools coverage -counts -a sliding_window.bed -b bed/KK_R_2.bed > count/KK_R_2.bed



