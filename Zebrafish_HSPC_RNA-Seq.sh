
thread=6
rawdata=/p300s/yangyg_group/sunbf/WangMK/SGG_sc-m5C/20240715_Berry_zberafish_HSPC_merge/zebrafish_HSPC_RNASeq
dir=/p300s/yangyg_group/wangmk/SGG_sc-m5C/RNA-Seq_zebrafish_HSC/20240715_RNA-Seq_merge

zv10=/xtdisk/yangyg_group/wangmk/reference/zebrafish/hisat2/danRer10
zv10_fa=/xtdisk/yangyg_group/wangmk/reference/zebrafish/zv10.fa
zv10_chr=/xtdisk/yangyg_group/wangmk/reference/zebrafish/ref/single-chr-trans
zv10_gtf=/xtdisk/yangyg_group/wangmk/reference/zebrafish/gtf_zv10/Danio_rerio.GRCz10.91.gtf

#fastqc
cd ${dir}
mkdir zebrafish_HSPC_RNASeq
cd ./zebrafish_HSPC_RNASeq

#��һ����clean data��
mkdir 1-clean_data
cd 1-clean_data

mkdir fastqc

fastqc $rawdata/*R1.fastq.gz -t $thread -o ./fastqc &
fastqc $rawdata/*R2.fastq.gz -t $thread -o ./fastqc &
wait

###ȥ��ͷ
cutadapt -a GATCGGAAGA -A GATCGGAAGA -o cutadapt_1.fq.gz -p cutadapt_2.fq.gz $rawdata/*R1.fastq.gz $rawdata/*R2.fastq.gz

fastqc cutadapt_1.fq.gz -t $thread -o ./fastqc &
fastqc cutadapt_2.fq.gz -t $thread -o ./fastqc &
wait

###����/�������
java -Xmx4g -jar /software/biosoft/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ./cutadapt_1.fq.gz ./cutadapt_2.fq.gz ./trim_1.fq.gz ./unpaired_1.fq.gz ./trim_2.fq.gz ./unpaired_2.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:35

rm unpaired_1.fq.gz unpaired_2.fq.gz

fastqc trim_1.fq.gz -t $thread -o ./fastqc &
fastqc trim_2.fq.gz -t $thread -o ./fastqc &
wait

#�ڶ����������rDNA
cd $dir/zebrafish_HSPC_RNASeq
mkdir 2-rDNA
cd 2-rDNA
rDNA=/xtdisk/yangyg_group/wangmk/reference/zebrafish/rDNA/zb_rDNA
hisat2 -p $thread --dta --rna-strandness RF -x $rDNA -1 ../1-clean_data/trim_1.fq.gz -2 ../1-clean_data/trim_2.fq.gz -S ./rDNA.sam --un-conc . > rDNA_mapping.txt 2>&1
rm rDNA.sam
mv un-conc-mate.1 de_rDNA_1.fastq
mv un-conc-mate.2 de_rDNA_2.fastq
gzip de_rDNA_*

#�ڶ��������бȶ�
####method1
cd $dir/zebrafish_HSPC_RNASeq
mkdir 3-genome
hisat2 -p $thread --dta --rna-strandness RF -x $zv10 -1 ./2-rDNA/de_rDNA_1.fastq.gz -2 ./2-rDNA/de_rDNA_2.fastq.gz -S ./3-genome/hisat2.sam --un-conc ./3-genome/ > ./3-genome/genome_mapping.txt 2>&1

rm ./3-genome/un-conc-mate*

####
cd 3-genome
###����ͳ��
#ͳ��δ�ʿص�reads��Ŀ
echo "reads with no -q:" > readsCount-genecode
cat hisat2.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode
#��bam�ļ������ʿ�
samtools view -hb hisat2.sam > hisat2.bam
samtools view hisat2.bam -q 20 -h > hisat2_20.sam
#ͳ���ʿغ��reads��Ŀ
echo "reads with -q:" >> readsCount-genecode
cat hisat2_20.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode
#��sam�ļ�����uniqueɸѡ
cat hisat2_20.sam|grep -E "^@|NH:i:1$|NH:i:1[^0-9]" > uniqmap.sam
#ͳ��unique��reads��Ŀ
echo "reads with uniqmap:" >> readsCount-genecode
cat uniqmap.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode
#ͳ��unique�������pair reads��Ŀ
echo "pair reads with uniqmap:" >> readsCount-genecode
cat uniqmap.sam|awk '$7~/=/'|cut -f1|awk '!x[$0]++'|wc -l >> readsCount-genecode
#תunique��sam�ļ�Ϊbam�ļ�
samtools view -S uniqmap.sam -b -o uniqmap.bam
#ɾ��sam�ļ�
rm hisat2_20.sam
rm hisat2.sam
rm uniqmap.sam

###����reads����ע�� 
#תbam�ļ�Ϊbed�ļ�
bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed
#���ݲο��ļ��õ�������Ϣ
intersectBed -a uniqmap.bed -b $zv10_chr -wa -wb -f 0.51 > tmp.SG2
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
#ͳ��ת¼��ıȶ�reads��
echo "InputRNA mapping reads:" >> readsCount-genecode
cat uniqmap.bed|awk '!x[$4]++'|wc -l >> readsCount-genecode
#ͳ�ƷǱ���RNA��reads��Ŀ
echo "InputRNA ncRNA reads:" >> readsCount-genecode
cat tmp.SG2 |awk '$14 != "protein_coding" && !a[$4]++'|wc -l >> readsCount-genecode
#ͳ��mRNA����Ŀ
echo "InputRNA mRNA reads:" >> readsCount-genecode
cat tmp.SG2 |awk '$14 == "protein_coding" && !a[$4]++'|wc -l >> readsCount-genecode
#ͳ������λ�������reads��Ŀ
echo "InputRNA all region reads:" >> readsCount-genecode
cat tmp.SG2|cut -f13|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount-genecode

echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
#ͳ�����бȶԵ���RNA���༰����
echo "InputRNA all RNA varieties:" >> readsCount-genecode
cat tmp.SG2|cut -f14|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount-genecode

echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode

#mRNA�ϵ�ע��
echo "InputRNA region reads only for mRNA:" >> readsCount-genecode
cat tmp.SG2|awk '$14 == "protein_coding"'|cut -f13|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode
echo "" >> readsCount-genecode

rm tmp.SG2
rm uniqmap.bed

#reads countͳ��
#����ÿ�������reads��������ѡhtSeq-count/StringTie�Ⱦ���
samtools sort uniqmap.bam -@ 4 > sort.bam
rm uniqmap.bam

featureCounts -T 4 -a $zv10_gtf -t exon -g gene_id -p -s 2 -Q 20 -R BAM -o InputRNA-featureCounts-s2-p.txt ./sort.bam
cat InputRNA-featureCounts-s2-p.txt|sed '1,2 d'|cut -f1,7|sort -k1 -n > featureCounts_filter.out

rm sort.bam.feature*.bam
#####
###����bw�ļ�
samtools index sort.bam

/software/biosoft/software/python/python2.7_2018_12/bin/bamCoverage --normalizeUsing RPKM -b sort.bam -o zebrafish_HSPC_RNASeq_sort.bw


