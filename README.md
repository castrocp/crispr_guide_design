# crispr_guide_design

Conda environment is set up to run the scripts here  
`conda activate chopchop`

## To run ChopChop script
options:  
-b (bedifle)  
-f (flank)  
-d (distance)  
-o (output)  
-t (threads)    

Ex.  
`sh nanopore_chopchop_bash.sh -b /home/castrocp/projects/guidedesign/48_targets.bed -d 1000 -f 200 -o ~/projects/guidedesign/TEST_output/ -t 2`  

### First part of script generates .txt files containing candidate guide sequences

### Second part filters the chopchop results
Criteria filtered on:\
min_gc and max_gc - GC content\
max_self_complementarity - number indicates how many regions of self-complementarity are predicted\
min_efficiency_score\
max_mm0/1/2/3 - hot many off-targets each target site has with 0 (mm0), 1 (mm1), 2 and 3 mismatches  

A .tsv is created for the filtered guides, and then a .bed file is created from that.  

### Third part creates fasta from guides, and scores them with crisprON

FASTA files are used as input for crisprON and output is a .csv with each guide scored.  

Guides are sorted by score from high to low and top 3 are sent to a `high_scoring_guides` BED file.

All high_scoring_guides are combined into a single `output_guides` BED file.
