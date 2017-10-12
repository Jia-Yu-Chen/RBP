#!/bin/bash -x

cell=${1:-HepG2}
window=${2:-2000000}
chrsize=${3:-~jchen/public/hg19.chrom.sizes}
cpg=${4:-~jchen/public/cpgIslandExt.hg19.txt}

# genome partition
more $chrsize | head -n 24 | bedtools makewindows -g /dev/stdin -w $window > genomeBins
# cpg island
bedtools intersect -a genomeBins -b $cpg -c > genomeBins.GC
maxgc=`sort -k 4,4nr genomeBins.GC | cut -f 4 | head -1`
awk -v max=$maxgc 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/max}' genomeBins.GC > a && mv a genomeBins.GC

# number of RBP peaks
ls ~jchen/RBP/data/$cell/narrowPeak/*narrowPeak ~jchen/RBP/data/$cell/narrowPeak/lowquality/*narrowPeak | while read f;
do
  cut -f 1-3 $f | sort -u
done | bedtools intersect -a genomeBins -b - -c | sed 's/^chr/hs/g'> genomeBins.RBPnum.${cell}

# Histone maker
if [ "$cell" = "HepG2" ]
then
  bigWigAverageOverBed ~jchen/RBP/data/HepG2/Epigenome/GSM816662_hg19_wgEncodeOpenChromDnaseHepg2Sig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > Dnase.HepG2
  bigWigAverageOverBed ~jchen/RBP/data/HepG2/Epigenome/wgEncodeBroadHistoneHepg2H3k4me3StdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K4me3.HepG2
  bigWigAverageOverBed ~jchen/RBP/data/HepG2/Epigenome/wgEncodeBroadHistoneHepg2H3k36me3StdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K36me3.HepG2
  bigWigAverageOverBed ~jchen/RBP/data/HepG2/Epigenome/wgEncodeBroadHistoneHepg2H3k27acStdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K27ac.HepG2
  bigWigAverageOverBed ~jchen/RBP/data/HepG2/Epigenome/wgEncodeBroadHistoneHepg2H3k09me3Sig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K9me3.HepG2
  cat H3K9me3.HepG2 | perl -ne '@t=split; $t[3]=log($t[3]) if $t[3]>0; $t[3] = 0 if $t[3] <=0; print join("\t",@t),"\n";' > a && mv a H3K9me3.HepG2
#  bigWigAverageOverBed ~jchen/RBP/data/HepG2/Epigenome/wgEncodeBroadHistoneHepg2CtcfStdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
#    cut -f 5 | paste genomeBins - | sed 's/^chr/hs/g'> Ctcf.HepG2
fi

if [ "$cell" = "K562" ]
then
  bigWigAverageOverBed ~jchen/RBP/data/K562/Epigenome/GSM816655_hg19_wgEncodeOpenChromDnaseK562SigV2.bigWig  <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > Dnase.K562
  bigWigAverageOverBed ~jchen/RBP/data/K562/Epigenome/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K4me3.K562
  bigWigAverageOverBed ~jchen/RBP/data/K562/Epigenome/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K36me3.K562
  bigWigAverageOverBed ~jchen/RBP/data/K562/Epigenome/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K27ac.K562
  bigWigAverageOverBed ~jchen/RBP/data/K562/Epigenome/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
    cut -f 5 | paste genomeBins - | sed 's/chr/hs/g' > H3K9me3.K562
  cat H3K9me3.K562 | perl -ne '@t=split; $t[3]=log($t[3]) if $t[3]>0; $t[3] = 0 if $t[3] <=0; print join("\t",@t),"\n";' > a && mv a H3K9me3.K562
#  bigWigAverageOverBed ~jchen/RBP/data/K562/Epigenome/wgEncodeBroadHistoneK562CtcfStdSig.bigWig <(awk 'BEGIN{OFS="\t"}{print $0,NR}' genomeBins) /dev/stdout | \
#    cut -f 5 | paste genomeBins - | sed 's/^chr/hs/g'> Ctcf.K562
fi
