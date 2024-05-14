#!/bin/bash
for f in /nfs/boylelab_turbo/nanopore_data/RepeatExpansion/20231102*/*.sorted.bam
do
echo "${f}"
#samtools stats ${f} | grep ^SN | cut -f 2- | grep 'bases mapped (cigar)'
samtools stats ${f} | grep ^SN | cut -f 2- | grep 'raw total sequences'
#samtools bedcov /home/crmumm/STR/STR_Panel/202302_targets.bed ${f} | cut -f 4
#samtools bedcov /home/crmumm/STR/guides/20210323_guides.sorted.bed ${f} | cut -f 4
samtools view -b -F 0x100 ${f} | bedtools coverage -a /home/crmumm/STR/STR_Panel/60_targets.bed -b stdin | cut -f 4
#bedtools coverage -a /home/crmumm/STR/guides/20210323_guides.sorted.bed -b ${f} | cut -f 4
done

