####macs2 call peak analysis

export PATH=/software/biosoft/software/python/python2.7_2018_12/bin:$PATH
dir1_1=/p300s/yangyg_group/wangmk/RIP-Seq_zebrafish_HSC/Input/3-genome
dir2_1=/p300s/yangyg_group/wangmk/RIP-Seq_zebrafish_HSC/RIP_rep1/3-genome
dir2_2=/p300s/yangyg_group/wangmk/RIP-Seq_zebrafish_HSC/RIP_rep2/3-genome

cd /p300s/yangyg_group/wangmk/SGG_sc-m5C/RIP-Seq_zebrafish_HSPC/macs2

mkdir macs2_callpeak
cd macs2_callpeak

mkdir keep-dup-p0.05
cd keep-dup-p0.05

macs2 callpeak -t $dir2_1/sort.bam $dir2_2/sort.bam -c $dir1_1/sort.bam --keep-dup all -f BAM --nomodel -g 1369631918 -B -p 0.05

cut -f1,2,3,7 NA_peaks.narrowPeak > ok.bed
perl /xtdisk/yangyg_group/wangmk/reference/zebrafish/ref/distribution_plot_zebrafish.pl

