pos=$(grep -oP '^MGAM,.*,7,\K\d*' SNPs_of_interest.csv)

#echo "$pos"

grep -f <(cat <<< "$pos") SNPs_of_interest_analysis_clst.frq.strat
