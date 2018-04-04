module load bedtools/2.21.0
out=/srv/gsfs0/projects/kundaje/users/oursu/other/junkDNA/mouse/mm9/
genefile=${out}genes
tssfile=${out}tss
cagefile=${out}cage
dnasefile=${out}dnase
histonefile=${out}histone
consfile=${out}cons
tffile=${out}tfs
humanfile=${out}human
repeatsfile=${out}mm9.repeatmasker.bed
blacklistfile=${out}mm9-blacklist.bed.gz
unmappable=${out}unmappable

#mm9 combine all coordinates for
#genes + predicted TSSs + CAGE peaks + dnase + histone + tfs + conserved

#### DNase
dnases=$(ls /srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/peaks_dnase/dcc/*/*)
for m in $dnases;
do
    echo $m
    zcat ${m} >> ${dnasefile}.tmp
done
cat ${dnasefile}.tmp | cut -f1-3 |  sort -k1,1 -k2,2n | bedtools merge -i - > ${dnasefile}
#### ==========================================================

### CONSERVED
#gerp
gerpdir=/srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/annotations/conservation/elements/mm9/GERP
for chromo in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y;
do
    echo $chromo
    zcat ${gerpdir}/mm9_chr${chromo}_elems.txt.gz | awk -v var="$chromo" '{print "chr"var"\t"$1"\t"$2}' >> ${out}gerp.tmp
done
#phastcons
phastcons=$(ls /srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/annotations/conservation/elements/mm9/phastCons/*)
for phasconsfile in ${phastcons};
do
    zcat ${phasconsfile} | cut -f2-4 >> ${out}phastcons
done
#phylo
phylos=$(ls /srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/annotations/conservation/elements/mm9/phyloP/*)
for phylo in ${phylos};
do
    zcat ${phylo} | cut -f2-4 >> ${out}phylo
done
cat ${out}gerp.tmp ${out}phastcons ${out}phylo | sort -k1,1 -k2,2n | bedtools merge -i - > ${consfile}
#####===============================================================

##### genes
gencode=/srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/annotations/gencodeAnno/mm9/vM1/gencode.vM1.annotation.gtf.gz
zcat ${gencode} | sed 's/;/\t/g' | awk '{if ($3=="UTR" || $3=="exon") print $1"\t"$4"\t"$5+1"\t"$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > ${genefile}
##### ===============================================================

##### tss
tss=/srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/annotations/tss/mm9/fantom5_tss/TSS_mouse.bed.gz
tss2=/srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/annotations/tss/mm9/gencode_tss/mm65_capped_sites_nr_with_confidence.gff.gz
zcat ${tss2} | sed 's/;/\t/g' | awk '{if ($3=="UTR" || $3=="exon") print $1"\t"$4"\t"$5+1"\t"$3}' | gzip > ${out}tss2.tmp.gz
zcat ${tss} ${out}tss2.tmp.gz | sort -k1,1 -k2,2n | bedtools merge -i - > ${tssfile}
#### ================================================================

#### cage
cagedir=/srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/peaks_CAGE/mm9/
zcat ${cagedir}mm9.cage_peak_coord_*.bed.gz | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > ${cagefile}
#### =================================================================

#### histone
histonedirfiles=$(ls /srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/peaks_histones/dcc/*/*)
for hfile in ${histonedirfiles};
do
    echo ${hfile}
    zcat ${hfile} >> ${out}histone.tmp
done
cat ${out}histone.tmp | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > ${histonefile}
#### =================================================================

#### TFs
tfdirfiles=$(ls /srv/gsfs0/projects/kundaje/commonRepository/mouseENCODE/peaks_spp/combrep/*)
for tfsppfile in ${tfdirfiles};
do
    echo ${tfsppfile}
    zcat ${tfsppfile} >> ${out}tf.tmp
done
cat ${out}tf.tmp | cut -f1-3 | awk '{printf $1"\t""%.0f\t%.0f\n",$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > ${tffile}
#### =========================================================================

#### unmappable
mapfile=${out}crgMapabilityAlign36mer.bw
/home/oursu/devtools/bigWigMerge -threshold 1.8 ${mapfile} ${mapfile} > ${mapfile}.mappable09.bedgraph
unmappable=${out}unmappable
mm9gen=/srv/gs1/software/bedtools/2.21.0/genomes/mouse.mm9.genome
cat ${mapfile}.mappable09.bedGraph | bedtools slop -b 0 -i - -g ${mm9gen} | bedtools complement -i - -g ${mm9gen} > ${unmappable}

### put together everything we have so far
cat ${tffile} ${histonefile} ${cagefile} ${tssfile} ${genefile} ${dnasefile} ${out}gerp.tmp ${out}phastcons | sort -k1,1 -k2,2n | bedtools merge -i -> ${out}_tf_hist_cage_tss_gene_dhs_cons

#Liftover all the combined human positive set to mouse mm9add to the above
#then take the complement
#and remove blacklist
human_set=${out}human_set___tss_gene_bl_unmap_dhs_repeat_cons_chipseq
hdir=/srv/gsfs0/projects/kundaje/users/oursu/other/junkDNA/
cat ${hdir}tss ${hdir}gencodeExonUTR ${hdir}unmappable ${hdir}dhs ${hdir}repeats ${hdir}conserved ${hdir}chip | sort -k1,1 -k2,2n | bedtools merge -i - > ${human_set}
/home/oursu/devtools/liftOver ${human_set} /home/oursu/devtools/liftOverData/hg19ToMm9.over.chain ${humanfile} ${humanfile}_unmapped

zcat ${blacklistfile} > ${blacklistfile}.tmp
cat ${repeatsfile} | cut -f1-3 > ${repeatsfile}.tmp
cat ${humanfile} ${repeatsfile}.tmp ${blacklistfile}.tmp ${unmappable} | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_human_repeats_blacklist_unmappable
cat ${out}_human_repeats_blacklist_unmappable ${out}_tf_hist_cage_tss_gene_dhs_cons | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_tf_hist_cage_tss_gene_dhs_cons_human_repeats_blacklist_unmappable


mm9gen=/srv/gs1/software/bedtools/2.21.0/genomes/mouse.mm9.genome
mm9_f=${out}_tf_hist_cage_tss_gene_dhs_cons_human_repeats_blacklist_unmappable
cat ${mm9_f} | bedtools slop -b 0 -i - -g ${mm9gen} | bedtools complement -i - -g ${mm9gen} | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_tf_hist_cage_tss_gene_dhs_cons_human_repeats_blacklist_unmappable_COMPLEMENT

cat ${out}_tf_hist_cage_tss_gene_dhs_cons_human_repeats_blacklist_unmappable_COMPLEMENT | awk '{if ($3-$2>2000) print $0}' | wc -l

#lift over these regions to mm10
mm9controls=${out}_tf_hist_cage_tss_gene_dhs_cons_human_repeats_blacklist_unmappable_COMPLEMENT
/home/oursu/devtools/liftOver ${mm9controls} /home/oursu/devtools/liftOverData/mm9ToMm10.over.chain ${mm9controls}liftOver2mm10.bed ${mm9controls}liftOver2mm10.unmapped 