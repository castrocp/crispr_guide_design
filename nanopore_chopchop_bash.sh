#!/bin/bash
source /home/castrocp/miniconda3/etc/profile.d/conda.sh
conda activate chopchop

##Changes with each experiment
while getopts b:f:d:o:t: flag
do
  case "$flag" in
    b) BEDFILE=${OPTARG};;
    f) FLANK=${OPTARG};;
    d) DISTANCE=${OPTARG};;
    o) OUTPUT=${OPTARG};;
    t) THREADS=${OPTARG};;
    *) echo "Invalid option: -$flag" exit 1;;
  esac
done

mkdir -p "$OUTPUT"

do_chopchop () {
  DISTANCE="$2"
  FLANK="$3"
  OUTPUT="$4"
  # # loop through each line in the file
  IFS=$'\t' read -r col1 col2 col3 <<< "$1"
  # make the second and third columns integer variables, if possible
  if [[ $col2 =~ ^[0-9]+$ ]]; then
    col2_int=$((col2))
  else
    col2_int="not an integer"
    echo 'stop coord "$col2" is "$col2_int"'
    exit
  fi
  col3="$(echo -e "${col3}" | sed -e 's/[[:space:]]*$//')"
  if [[ $col3 =~ ^[0-9]+$ ]]; then
    col3_int=$((col3))
  else
    col3_int="not an integer"
    echo 'stop coord "$col3" is "$col3_int'
    exit
  fi
  ############################## Makes regions and pass to chopchop

  upstream_coods_start=$(($col2_int - $DISTANCE - $FLANK))
  upstream_coods_end=$(($col2_int - $DISTANCE + $FLANK))

  downstream_coods_start=$(($col3_int + $DISTANCE-$FLANK))
  downstream_coods_end=$(($col3_int + $DISTANCE+$FLANK))

  LOC="$col1"."$col2_int"."$col3_int"

  # echo "$col1":"$upstream_coods_start"-"$upstream_coods_end"
  # echo "$col1":"$downstream_coods_start"-"$downstream_coods_end"

  #Upstream
  python /home/castrocp/software/chopchop/chopchop.py \
    -Target "$col1":"$upstream_coods_start"-"$upstream_coods_end" \
    -BED \
    -GenBank \
    -G hg38 \
    -filterGCmin 20 \
    -filterGCmax 80 \
    -filterSelfCompMax 0 \
    -t WHOLE \
    -n N \
    -R 4 \
    -T 1 \
    -g 20 \
    -scoringMethod DOENCH_2016 \
    -f NN \
    -v 3 \
    -M NGG \
    -BB AGGCTAGTCCGT \
    -o "$OUTPUT"/"$LOC"_temp \
    --nonOverlapping | awk '$4 == "+" {print}' | head -n 25 >> "$OUTPUT"/"$LOC".upstream.txt

  #Downstream
  python /home/castrocp/software/chopchop/chopchop.py \
    -Target "$col1":"$downstream_coods_start"-"$downstream_coods_end" \
    -BED \
    -GenBank \
    -G hg38 \
    -filterGCmin 20 \
    -filterGCmax 80 \
    -filterSelfCompMax 0 \
    -t WHOLE \
    -n N \
    -R 4 \
    -T 1 \
    -g 20 \
    -scoringMethod DOENCH_2016 \
    -f NN \
    -v 3 \
    -M NGG \
    -BB AGGCTAGTCCGT \
    -o "$OUTPUT"/"$LOC"_temp \
    --nonOverlapping | awk '$4 == "-" {print}' | head -n 25 >> "$OUTPUT"/"$LOC".downstream.txt
  rm -r "$OUTPUT"/"$LOC"_temp
}

export -f do_chopchop

parallel -j "$THREADS" -d "\n" do_chopchop {} "$DISTANCE" "$FLANK" "$OUTPUT" :::: "$BEDFILE"

############################### Make 30mers

#for candidate in "$OUTPUT"/*.txt
#do
#
#name=$(basename $candidate .txt)
#
## Filter criteria
#min_gc=20
#max_gc=80
#max_self_complementarity=0
#min_efficiency_score=0.3
#max_mm0=0
#max_mm1=1
#max_mm2=1
#max_mm3=10
#
## Read input file line by line and filter based on criteria
#while IFS=$'\t' read -r id sequence location strand gc selfcomp mm0 mm1 mm2 mm3 efficiency_score; do
#    if (( gc >= min_gc && gc <= max_gc )); then
#      if ((selfcomp <= max_self_complementarity)); then
#        if [[ $mm0 != *">="* && $mm1 != *">="* && $mm2 != *">="* && $mm3 != *">="*  ]]; then
#            if (( mm0 <= max_mm0 && mm1 <= max_mm1 && mm2 <= max_mm2 ** mm3 <= max_mm3)); then
#              if (( $(echo "$efficiency_score > $min_efficiency_score" | bc -l) )); then
#                echo -e "$id\t$sequence\t$location\t$strand\t$gc\t$selfcomp\t$mm0\t$mm1\t$mm2\t$mm3\t$efficiency_score" >> "$OUTPUT"/"$name"_filtered.tsv
#              fi
#            fi
#          fi
#        fi
#      fi
#done < "$candidate"
#
#done
#
#for candidate in "$OUTPUT"/*_filtered.tsv
#do
#
#name=$(basename $candidate _filtered.tsv)
#
#while read line; do
#
#chr=$(echo "$line" | awk '{split($3,a,":"); print a[1]}')
#position=$(echo "$line" | awk '{split($3,a,":"); print a[2]}')
#strand=$(echo "$line" | awk '{print $4}')
## Build bed for guides
#
#if [[ "$strand" = "+" ]]
#  then
#    mer_start=$(("$position" - 5))
#    mer_end=$(("$position" + 25))
#    echo "$chr" "$mer_start" "$mer_end" "." "." "$strand" | sed 's/ /\t/g' >> "$OUTPUT"/"$name".bed
#    
#  else
#    mer_start=$(("$position" - 4))
#    mer_end=$(("$position" + 26))
#    echo "$chr" "$mer_start" "$mer_end" "." "." "$strand" | sed 's/ /\t/g' >> "$OUTPUT"/"$name".bed
#
#fi
#
##echo "$chr" "$mer_start" "$mer_end" "." "." "$strand" | sed 's/ /\t/g' >> "$OUTPUT"/"$name".bed
#
## Done inner while loop
#done < "$candidate"
#
## Done
#done

# #################################### Make 30mer fastas and run


# WILL STILL WANT TO RUN EVERYTHING BELOW THIS, BUT FOR NOW JUST WANT TO RUN CHOPCHOP

#conda deactivate
#conda activate crispron
#
#do_crispron() {
#  candidate_bed=$1
#  OUTPUT=$2
#  POS=$3
#  name=$(basename $candidate_bed .bed)
#
#  # Use bed coords to make fastas for scoring
#  /home/crmumm/software/bedtools getfasta -s -fi /home/crmumm/data/hg38.fa -bed $candidate_bed > "$OUTPUT"/"$name"."$POS".fasta
#
#  #Score with crisprON TL
#  /home/xinleee/TL-crispron/bin/CRISPRon.sh "$OUTPUT"/"$name"."$POS".fasta "$OUTPUT"/"$name"
#
#  mv "$OUTPUT"/"$name"/crispron.csv "$OUTPUT"/"$name"_crispron_results."$POS".csv
#  rm -r "$OUTPUT"/"$name"
#  # Gets the 3 highest scoring guides - no filter for score > 10 
#  sort -k3 -nr -t, "$OUTPUT"/"$name"_crispron_results."$POS".csv > "$OUTPUT"/"$name"_crispron_results."$POS".sorted.csv
#  head -n -1 "$OUTPUT"/"$name"_crispron_results."$POS".sorted.csv | awk -F'[(),]|[_]' 'NR < 4 {split($1,a,":|-"); print a[1]"\t"a[2]"\t"a[3]"\t"$6"\t"$7"\t"$2}' > "$OUTPUT"/"$name"_high_scoring_guides."$POS".bed
#}
#
#export -f do_crispron
#
#ls "$OUTPUT"/*.upstream.bed| parallel -j "$THREADS" do_crispron {} "$OUTPUT" upstream
#ls "$OUTPUT"/*.downstream.bed | parallel -j "$THREADS" do_crispron {} "$OUTPUT" downstream
#
#cat "$OUTPUT"/*high_scoring_guides.*bed > "$OUTPUT"/output_guides.bed

