# Various-script-2023-master-thesis
Various scripts and results I produced during my master thesis at the Center of Evolutionnary Biology, in Uppsala, Sweden.

Details of each file will be added later

## Retrieve genes

The first scripts that were made were those for retrieving the genes of interest. By retrieving, I mean multiple things :  

For the differents SNPs positions that are presents in the H3Africa SNP chip :
- Filter positions of the SNPs associated to the genes of interest in the list of SNPs presents in the H3Africa micro array. This was done with the script [retrieve_gene.ps1](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/retrieve_genes.ps1) which is just a grep command looking for any line containing the name of our genes of interest.
- Retrieve the rsID and alleles of each SNPs we got with the previous step. This was done with the script [dbSNP_retrieve.py](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/dbSNP_retrieve.py). This script opens a csv file containing the previous SNPs, and makes queries on NCBI Entrez SNP database to get the information that we want for each SNP. Then, it writes the results in a "dbSNP_retrieved.csv" CSV file. Doesn't take any parameters but the code opens a file with a specific name.

For the actual data genotyped from different populations in Africa, Europe and EastAsia :
- Retrieve every SNP of interest by using the positions of the SNPs that we got with the step #1. This was done with the scripts [get_positions.sh](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/grep_positions.sh) and [sbatch_grep_positions.sh](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/sbatch_grep_positions.sh), which are just a grep command looking for lines containing the positions of our SNPs.

## Caracterizing genetic diversity

Thanks to the [plink](https://www.cog-genomics.org/plink/) software, I computed the frequency of each SNP in each population with the command

```bash
./plink --freq --family --tfile SNPs_of_interest --out SNPs_of_interest_analysis_clst
```

### [modify_duplicate_coords.py](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/modify_duplicate_coords.py)

This basic python scripts checks in the file with the populations' coordinates if there are duplicates coordinates, and if in this case slightly changes the coords.

### [pie_charts_allele_frequencies.R](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/pie_charts_allele_frequencies.R)

This R script takes the file containing populations' coordinates (corrected with `modifiy_duplicate_coords.py`) and the file with SNPs frequencies, and represent them on a world map. An exemple of map is visible at the adress : https://paul-bunel.github.io/Various-script-2023-master-thesis/map.html

### [analyzes.R](https://github.com/Paul-bunel/Various-script-2023-master-thesis/blob/main/analyzes.R)

This R script takes the same files as the previous one, and contains differents commands to make a table of each SNPs frequency for each population. Then different basics plots are implemented, in particular a plot showing the mean frequency of each SNP for the three main population's lifestyles categories : farmers, pastoralists and foragers. This plot is visible at the adress : https://paul-bunel.github.io/Various-script-2023-master-thesis/plot.html

## Testing selection

Thanks to the [plink2](https://www.cog-genomics.org/plink/2.0/) software, I computed the pairwise $F_{ST}$ statistic for some groups of population with the command

```bash
../plink2 --fst CATPHENO report-variants --within clst_pops_of_interest.tfam --tfile ../AfricanNeo_and_Public_DB_minN10_Jan2022 --keep-fam pops_of_interest.tfam --out pops_of_interest_fst
```

This command use the files `pops_of_interest.tfam`, that contains a subset of the main tfam file with only the population of interest, and `clst_pops_of_interest.tfam`, that is the exact same file with cluster names added at the end of each sample line. These files were generated with the following bash commands

```bash
grep -wP '.*Japanese.*|.*Baka.*|.*Nigeria_Yoruba.*|DRC_((?!LubaLulua|Shi|Rega).)*' AfricanNeo_and_Public_DB_minN10_Jan2022.tfam > pops_of_interest/pops_of_interest.tfam

sed -r "s/(DRC.*\s.*)\s\w\s\w\s\w\s\w/\1 DRC/" pops_of_interest.tfam > tmp.tfam
sed -r "s/(.*Baka\s.*)\s\w\s\w\s\w\s\w/\1 Baka/" tmp.tfam > clst_pops_of_interest.tfam
sed -r "s/(.*Yoruba.*\s.*)\s\w\s\w\s\w\s\w/\1 Yoruba/" clst_pops_of_interest.tfam > tmp.tfam
sed -r "s/(.*Japanese.*\s.*)\s\w\s\w\s\w\s\w/\1 Japanese/" tmp.tfam > clst_pops_of_interest.tfam
```
