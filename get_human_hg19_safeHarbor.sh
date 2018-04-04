
######### Get the quiescent regions to start from ===================================================================================================================
module load bedtools/2.21.0
out=/srv/gsfs0/projects/kundaje/users/oursu/other/junkDNA/
quies=${out}quiescentIntersect
consfile=${out}conserved
dhsfile=${out}dhs
repeatsfile=${out}repeats
tffile=${out}TFs
blacklistfile=${out}blacklist
unmappablefile=${out}unmappable
gencodeExonUTRfile=${out}gencodeExonUTR
tssfile=${out}tss
chipfile=${out}chip

#### GET QUIESCENT STATE across all cell types in Roadmap 
roadmap=/srv/gsfs0/projects/kundaje/commonRepository/epigenomeRoadmap/segmentations/models/coreMarks/parallel/set2/final
mnemonics_files=$(ls ${roadmap}/*_15_coreMarks_mnemonics.bed.gz)
zcat ${roadmap}/E129_15_coreMarks_mnemonics.bed.gz | grep "Quies\|ReprPC\|ReprPCWk\|Het" > ${out}quiescentIntersect
cd ${roadmap}
for m in $mnemonics_files;
do
    echo $m
    zcat ${m} | grep "Quies\|ReprPC\|ReprPCWk\|Het" | bedtools intersect -a ${out}quiescentIntersect -b - > ${out}quiescentIntersect.tmp
    mv ${out}quiescentIntersect.tmp ${out}quiescentIntersect
    wc -l ${out}quiescentIntersect
done
cat ${out}quiescentIntersect | gzip > ${out}quiescentIntersect.gz
#### ============================================


#### Combine all conserved elements
consDir=/srv/gsfs0/projects/kundaje/commonRepository/annotations/human/hg19.GRCh37/conservedElements/
conservedElements=$(ls ${consDir}/*conserved.bed.gz)
cd ${consDir}
for m in $conservedElements;
do
    echo $m
    zcat ${m} >> ${out}conserved
done
cat ${out}conserved| sort -k1,1 -k2,2n | bedtools merge -i - > ${out}conserved.merged
mv ${out}conserved.merged ${out}conserved

#combine conserved domains too
conservedElements=$(ls ${consDir}/*domain.bed.gz)
cd ${consDir}
for m in $conservedElements;
do
    echo $m
    zcat ${m} >> ${out}conserved_domain
done
cat ${out}conserved_domain | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}conserved_domain.merged
mv ${out}conserved_domain.merged ${out}conserved_domain

#### ======================================


### Combine all DHSs
dhsdir=/srv/gsfs0/projects/kundaje/users/akundaje/projects/encode/data/byDataType/peaks_dnase/Uw/DNaseI_peaks/
dhsThings=$(ls ${dhsdir}*bed.gz)
cd ${dhsdir}
for m in $dhsThings;
do
    echo $m
    zcat ${m} >> ${dhsfile}.tmp
done
cat ${dhsfile}.tmp | sort -k1,1 -k2,2n | bedtools merge -i - > ${dhsfile}
#### ========================================


###### REPEATS
repdir=/srv/gsfs0/projects/kundaje/commonRepository/annotations/human/hg19.GRCh37/repeats/hg19/data/bed/
repfiles=$(ls ${repdir})
cp ${curfile} ${curfile}noRepeats
cd ${repdir}
for m in $repfiles;
do
    echo $m
    zcat ${m} >> ${repeatsfile}.tmp
done
cat ${repeatsfile}.tmp | sort -k1,1 -k2,2n | bedtools merge -i - > ${repeatsfile}
#### ========================================

###### TFBS. this is very large, so i will do some intermediate mergebed on the way
TFBS_encode=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/bindingSites/encodeMotifs/dec2013/globalHitLits/matches.txt.gz
zcat ${TFBS_encode} | split -l 1000000 - ${out}tf_chunk
tf_chunks=$(ls ${out}tf_chunk*)
for tf_chunk in ${tf_chunks};
do
    echo ${tf_chunk}
    cat ${tf_chunk} | sed 's/ /\t/g' | awk '{print $2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -i - > ${tffile}_$(basename ${tf_chunk})
done
cat ${tffile}_* | sort -k1,1 -k2,2n | bedtools merge -i - > ${tffile}
#### ========================================

##### CHIP-SEQ
chipdir=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/peaks_spp/mar2012/distinct/idrOptimalBlackListFilt/
chips=$(ls ${chipdir}*Tfbs*)
cd ${chipdir}
for chip in $chips;
do
    echo $chip
    zcat $chip >> ${chipfile}.tmp
done
cat ${chipfile}.tmp | sort -k1,1 -k2,2n | bedtools merge -i - > ${chipfile}

##### BLACKLIST and unmappable
blacklist=/srv/gs1/projects/kundaje/oursu/Alignment/data/blacklist/wgEncodeDacMapabilityConsensusExcludable.bed
zcat ${blacklist}.gz > ${blacklistfile}
unmappable=/srv/gsfs0/projects/kundaje/commonRepository/mappability/Park_lab/unmappable.hg19.anshul.25.bed.gz
zcat ${unmappable} > ${unmappablefile}
#### =========================================

#### Genes 
gencode=/srv/gs1/projects/snyder/jzaugg/histoneQTL/RNAseq/data/GENCODE_v19_2014-06-03/gencode.v19.annotation.gtf.gz
zcat ${gencode} | sed 's/;/\t/g' | awk '{if ($3=="UTR" || $3=="exon") print $1"\t"$4"\t"$5"\t"$3}' | sort -k1,1 -k2,2n | bedtools merge -i - > ${gencodeExonUTRfile}
#### ===========================================

###### TSSs                                                                                                                                      
TSS1=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/annotations/gencodeTSS/v19/fantom5.cage.permissive.hg19.bed.gz
TSS2=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/annotations/gencodeTSS/v19/fantom5.cage.strict.hg19.bed.gz
TSS3=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/annotations/gencodeTSS/v19/gencode.v19.annotation_capped_sites_nr_with_confidence.gff.gz
TSS4=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/annotations/gencodeTSS/v19/TSS_human_with_gencodetss_notlow_ext50eachside_merged_withgenctsscoord_andgnlist.gff.gz
TSS5=/srv/gsfs0/projects/kundaje/commonRepository/encode/data/byDataType/annotations/gencodeTSS/v19/TSS_human_strict_with_gencodetss_notlow_ext50eachside_merged_withgenctsscoord_andgnlist.gff.gz
for tss in $TSS1 $TSS2;
do
    echo $tss
    zcat ${tss} | sed 's/ /\t/g' >> ${tssfile}.tmp
done
cat ${tssfile}.tmp | grep -v "track" | sort -k1,1 -k2,2n | cut -f1-3 |  bedtools merge -i - > ${tssfile}
###### ========================================

#Now, download everything, and check where we are losing most regions
module load bedtools/2.21.0
out=/srv/gsfs0/projects/kundaje/users/oursu/other/junkDNA/
quies=${out}quiescentIntersect
consfile=${out}conserved
dhsfile=${out}dhs
repeatsfile=${out}repeats
tffile=${out}TFs
blacklistfile=${out}blacklist
unmappablefile=${out}unmappable
gencodeExonUTRfile=${out}gencodeExonUTR
tssfile=${out}tss

cat ${quies} | sort -k1,1 -k2,2n | bedtools merge -i - > ${quies}.merged

#remove tss, blacklist, unmappable, dhs, repeats
cat $tssfile $blacklistfile $unmappablefile $dhsfile $repeatsfile $gencodeExonUTRfile | sort -k1,1 -k2,2n | cut -f1-3 |  bedtools merge -i - > ${out}tss_gene_bl_unmap_dhs_repeat
bedtools subtract -a ${quies}.merged -b ${out}tss_gene_bl_unmap_dhs_repeat > ${out}_no_tss_gene_bl_unmap_dhs_repeat
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat | awk '{if ($3-$2>2000) print $0}' | wc -l
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat | awk '{if ($3-$2>500) print $0}'| wc -l

#remove conserved
bedtools subtract -a ${out}_no_tss_gene_bl_unmap_dhs_repeat -b ${consfile} > ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons | awk '{if ($3-$2>2000) print $0}'| wc -l
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons| awk '{if ($3-$2>500) print $0}'| wc -l

#remove tfbs motfis
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons.merged
bedtools subtract -a ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons -b ${tffile} > ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_tf
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_tf | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_tf.merged
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_tf.merged | awk '{if ($3-$2>2000) print $0}'| wc -l
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_tf.merged | awk '{if ($3-$2>500) print $0}'| wc -l

#instead of motifs, just remove chip-seq
bedtools subtract -a ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons -b ${chipfile} | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq | awk '{if ($3-$2>2000) print $0}'| wc -l
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq| awk '{if ($3-$2>500) print $0}'| wc -l

bedtools subtract -a ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq -b ${out}conserved_domain | sort -k1,1 -k2,2n | bedtools merge -i - > ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq_consdomain
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq_consdomain | awk '{if ($3-$2>2000) print $0}'| wc -l
cat ${out}_no_tss_gene_bl_unmap_dhs_repeat_cons_chipseq_consdomain | awk '{if ($3-$2>500) print $0}'| wc -l
