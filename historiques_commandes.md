Historiques commandes
=====================

Plink
-----

### Retrieve allelic frequencies by family

```bash
./plink --freq --family --tfile SNPs_of_interest --out SNPs_of_interest_analysis_clst
```

### Get population of interest

```bash
grep -wP '.*Japanese.*|.*Baka.*|.*Nigeria_Yoruba.*|DRC_((?!LubaLulua|Shi|Rega).)*' AfricanNeo_and_Public_DB_minN10_Jan2022.tfam > pops_of_interest/pops_of_interest.tfam

sed -r "s/(DRC.*\s.*)\s\w\s\w\s\w\s\w/\1 DRC/" pops_of_interest.tfam > tmp.tfam
sed -r "s/(.*Baka\s.*)\s\w\s\w\s\w\s\w/\1 Baka/" tmp.tfam > clst_pops_of_interest.tfam
sed -r "s/(.*Yoruba.*\s.*)\s\w\s\w\s\w\s\w/\1 Yoruba/" clst_pops_of_interest.tfam > tmp.tfam
sed -r "s/(.*Japanese.*\s.*)\s\w\s\w\s\w\s\w/\1 Japanese/" tmp.tfam > clst_pops_of_interest.tfam
```

### ? Fst computation ?

```bash
../plink2 --fst CATPHENO report-variants --within clst_pops_of_interest.tfam --tfile ../AfricanNeo_and_Public_DB_minN10_Jan2022 --keep-fam pops_of_interest.tfam --out pops_of_interest_fst
```