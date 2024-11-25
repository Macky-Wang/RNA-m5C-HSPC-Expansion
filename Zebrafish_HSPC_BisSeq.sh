thread=6
rawdata=/p300s/yangyg_group/sunbf/WangMK/SGG_sc-m5C/20230710_Berry_merge/zebrafish_HSPC_BisSeq
path=/p300s/yangyg_group/wangmk/SGG_sc-m5C/RNA-BisSeq_zebrafish_HSC/20230710_BisSeq_HSC_merge

zv10_fa=/xtdisk/yangyg_group/wangmk/reference/zebrafish/zv10.fa
zv10_chr=/xtdisk/yangyg_group/wangmk/reference/zebrafish/ref/single-chr-trans
zv10_gtf=/xtdisk/yangyg_group/wangmk/reference/zebrafish/gtf_zv10/Danio_rerio.GRCz10.91.gtf

cd $path
mkdir zebrafish_HSPC_BisSeq
cd zebrafish_HSPC_BisSeq
#####clean data
mkdir 1-clean_data
cd 1-clean_data
mkdir fastqc

##原始数据QC，查看质量
fastqc $rawdata/*R1.fastq.gz -t $thread -o ./fastqc &
fastqc $rawdata/*R2.fastq.gz -t $thread -o ./fastqc &
wait

cutadapt -a AGATCGGAAG -A AGATCGGAAG -o cutadapt_1.fastq.gz -p cutadapt_2.fastq.gz $rawdata/*R1.fastq.gz $rawdata/*R2.fastq.gz

java -Xmx4g -jar /pnas/yangyg_group/zhangting/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 cutadapt_1.fastq.gz cutadapt_2.fastq.gz ./trim_qua_1.fastq.gz ./unpaired_1.fastq.gz ./trim_qua_2.fastq.gz ./unpaired_2.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35
rm unpaired_*

fastqc ./trim_qua_1.fastq.gz -t $thread -o ./fastqc &
fastqc ./trim_qua_2.fastq.gz -t $thread -o ./fastqc &
wait

gunzip -c trim_qua_1.fastq.gz > trim_qua_1.fastq
gunzip -c trim_qua_2.fastq.gz > trim_qua_2.fastq

###rDNA filter
###zebrafish rDNA mapping
cd ${path}/zebrafish_HSPC_BisSeq
mkdir 2-rDNA_check
cd 2-rDNA_check
rDNA_meRanGh=/xtdisk/yangyg_group/wangmk/reference/zebrafish/rDNA/meRanGh_index
meRanGh align -o ./meRanGhResult -f $path/zebrafish_HSPC_BisSeq/1-clean_data/trim_qua_2.fastq -r $path/zebrafish_HSPC_BisSeq/1-clean_data/trim_qua_1.fastq -t ${thread} -S tmp.sam -un -ds -id $rDNA_meRanGh -mbp -fmo -mmr 0.01 > rDNA_meRanGh_ds.log 2>&1
rm meRanGhResult/tmp.bam
gzip meRanGhResult/tmp.sam

###rRNA filter
###zebrafish 45s RNA mapping
cd ${path}/zebrafish_HSPC_BisSeq
mkdir 3-rRNA_check
cd 3-rRNA_check
zv10_45s_meRanGh="/pnas/yangyg_group/yangxin/reference/zebrafish/45S/hisat2_index"
meRanGh align -o ./meRanGhResult -f $path/zebrafish_HSPC_BisSeq/2-rDNA_check/meRanGhResult/trim_qua_2_unmapped.fastq -r $path/zebrafish_HSPC_BisSeq/2-rDNA_check/meRanGhResult/trim_qua_1_unmapped.fastq -t ${thread} -S tmp.sam -un —ds -id $zv10_45s_meRanGh -mbp -fmo -mmr 0.01 > rRNA_meRanTK.log 2>&1
rm meRanGhResult/tmp.sam
fastqc trim_qua_1_unmapped_unmapped.fastq -t 6 -o ../../1-clean_data/fastqc/.
fastqc trim_qua_2_unmapped_unmapped.fastq -t 6 -o ../../1-clean_data/fastqc/.

gzip $path/zebrafish_HSPC_BisSeq/2-rDNA_check/meRanGhResult/trim_qua_*_unmapped.fastq

###genome mapping
cd $path/zebrafish_HSPC_BisSeq
mkdir 4-meRanGh
cd 4-meRanGh
meRanTK_index=/xtdisk/yangyg_group/wangmk/reference/zebrafish/meRanGh_index

meRanGh align -o ./meRanGhResult -f $path/zebrafish_HSPC_BisSeq/3-rRNA_check/meRanGhResult/trim_qua_2_unmapped_unmapped.fastq -r $path/zebrafish_HSPC_BisSeq/3-rRNA_check/meRanGhResult/trim_qua_1_unmapped_unmapped.fastq -t $thread -S tmp.sam -un -ds -id $meRanTK_index -GTF $zv10_gtf -mbp -fmo -mmr > genome_meRanGh.log 2>&1

cd meRanGhResult
rm *_unmapped.fastq
 
#####m5C calling
meRanCall -f $zv10_fa -p $thread -bam  tmp_sorted.bam -gref -o  m5C_call.txt -mBQ 20 -mr 0 >genome_calling.txt 2>&1

##m5C calling 10-2
cd ..
mkdir m5C_Calling_10-2
cd m5C_Calling_10-2

##考虑到mutation碱基，从总的覆盖碱基数目中去除突变碱基数目
cat ../meRanGhResult/m5C_call.txt | awk -v OFS="\t" '{if($12=="M")print $1,$2,$2,$3,$6,$5,$7,$12}' > methyC.bed

##注释并提取甲基化位点，合并后是覆盖该位置的总碱基数目（去除了mutation碱基的数目）
##这部分的注释，请根据研究课题来进行，不是一成不变的
cat methyC.bed|awk -v OFS="\t" '$7>=0.1 && $6>=10 && $5>=2' > methyC_qua.bed
intersectBed -a methyC_qua.bed -b $zv10_chr -wa -wb | awk '$4==$14'|awk '!a[$0]++' > methyC_gene.txt
cat methyC_gene.txt|awk '$16=="protein_coding"' > methy_C_mRNA.txt
cat methy_C_mRNA.txt|sed 's/_/\t/'|awk -v OFS="\t" '{print $1,$2,$3,$7,$15,$4}'|awk '!a[$1"\t"$2]++' > methylated.bed
cat methy_C_mRNA.txt|cut -f12|sort|uniq > 5mC_mRNA_genelist.txt
awk '{print $16}' methyC_gene.txt | sort | uniq -c | sort -k1 -nr | awk -v OFS="\t" '{print $2,$1}'> geneType.txt

sortBed -i methylated.bed | bedtools merge -i - -c 6 -o collapse -s -d 10 | grep -v "," | awk 'NR==FNR{a[$1_$2]=$0}NR>FNR{if(a[$1_$2]){print a[$1_$2]}}' methylated.bed - > filter_methylated.bed
mkdir filter_mRNA_distribution
cut -f1-3,6 filter_methylated.bed > filter_mRNA_distribution/distribution_mRNA.bed
cd filter_mRNA_distribution/
perl /xtdisk/yangyg_group/wangmk/reference/zebrafish/ref/distribution_plot_zv10.pl
cd ..
awk 'NR==FNR{a[$1_$2]=$0}NR>FNR{if(a[$1_$2]){print a[$1_$2]}}' methyC_gene.txt filter_methylated.bed | cut -f12 | awk '!a[$0]++' > filter_5mC_mRNA_genelist.txt

mkdir motif
cd motif

cat ../filter_methylated.bed | awk -v OFS="\t" '{print $1,$2-11,$3+10,"*","*",$6}'|fastaFromBed -fi $zv10_fa -bed - -s -fo ./motif-20.fasta
cat motif-20.fasta |grep -v ">"|sed 's/T/U/g' > weblogo-20.fa

###Luciferase mapping
cd $path/zebrafish_HSPC_BisSeq
mkdir 5-lucif
cd 5-lucif
htseq_index=/pnas/yangyg_group/yangxin/reference/Luciferase/Luciferase.fa
meRanGh_Lucif_index=/pnas/yangyg_group/yangxin/reference/Luciferase/meRanGh1.2.1_Luciferase

meRanGh align -o ./meRanGhResult -f ../3-rRNA_check/meRanGhResult/trim_qua_2_unmapped_unmapped.fastq -r ../3-rRNA_check/meRanGhResult/trim_qua_1_unmapped_unmapped.fastq -t 6 -S ./meRanGhResult/align.sam -un -ds -id $meRanGh_Lucif_index -mbp -fmo -mmr 0.01 > Lucif_mapping.txt 2>&1 
rm ./meRanGhResult/trim_qua_1_unmapped_unmapped_unmapped.fastq
rm ./meRanGhResult/trim_qua_2_unmapped_unmapped_unmapped.fastq

meRanCall -f $htseq_index -p 8 -bam ./meRanGhResult/align_sorted.bam -gref -o ./meRanGhResult/call.txt -mBQ 20 -mr 0 >Lucif_calling.txt 2>&1
rm ./meRanGhResult/align_sorted.bam*

gzip $path/zebrafish_HSPC_BisSeq/3-rRNA_check/meRanGhResult/trim_qua_*_unmapped_unmapped.fastq


rm $path/zebrafish_HSPC_BisSeq/1-clean_data/trim_qua_1.fastq
rm $path/zebrafish_HSPC_BisSeq/1-clean_data/trim_qua_2.fastq


