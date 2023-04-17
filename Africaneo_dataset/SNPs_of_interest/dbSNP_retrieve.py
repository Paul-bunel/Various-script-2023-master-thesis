from Bio import Entrez
from Bio.Data import IUPACData
import csv
from bs4 import BeautifulSoup
import time

# TODO : improve allele data retrieving by getting ancestral et derived alleles

# Reading CSV file with SNP position and Name, and saving it into dict with the
# following syntax : {position: {'gene': gene, 'ch': chromosome}}
snps_dict = {}
with open("retrieved_genes.csv", newline="") as csvfile:
    r = csv.reader(csvfile)
    for row in r:
        snps_dict[row[2].strip()] =\
            {"gene": row[0].strip(), "ch": row[1].strip()}

hg_version = "GRCh37"
API_KEY = "1b98794861449e649a073db4f086aa8b6308"

Entrez.email = "paul.bunel@ebc.uu.se"

# Entrez query to retrieve all NCBI UIDs of the SNPs in "retrieved_genes.csv"
# Query use position of the SNP, and the result (the UID) is added to the previous
# dict, the new syntax being : {position: {'gene': gene, 'ch': chromosome, 'uid': uid}}
for k, v in snps_dict.items():
    query = k + "[POSITION_GRCH37] AND " + v['ch'] + "[CHR]"
    handle = Entrez.esearch(db="snp", term=query, api_key=API_KEY)
    record = Entrez.read(handle)
    handle.close()
    if int(record['Count']) > 0:
        snps_dict[k]['uid'] = record['IdList'][-1]
    else:
        snps_dict[k]['uid'] = 0

uids = [snps_dict[k]['uid'] for k in snps_dict.keys()]
print("Len uids =", len(uids))

handle = Entrez.efetch(db='snp', id=uids, rettype="xml", retmode="xml")
records = handle.read()
soups = BeautifulSoup(records, "xml")
handle.close()

result_file = open("dbSNP_retrieved.csv", "w", newline="")
w = csv.writer(result_file)
for k, v in snps_dict.items():
    if v['uid'] != 0:
        s = soups.find_all(attrs={"uid": v['uid']})
        if len(s) > 1:
            print("=== ATTENTION PLUSIEURS OCCURENCES DE L'UID '%s'" % v['uid'])
        soup = s[0]

        try:
            gene = soup.NAME.get_text()
        except Exception as e:
            gene = "H3Africa : " + v['gene']
            print("gene name error :", e)
        rsID = "rs"+v['uid']
        try:
            ch, pos = soup.CHRPOS_PREV_ASSM.get_text().split(":")
        except Exception as e:
            ch = v['ch']
            pos = k
            print("ch/pos error:", e)
        try:
            allele_IUPAC = soup.ALLELE.get_text()
            decoded_ALLELE =\
                '|'.join([a for a in IUPACData.ambiguous_dna_values[allele_IUPAC]])
        except Exception as e:
            decoded_ALLELE = "N/A"
            print("allele error :", e)
        if gene != v['gene']  and gene != "H3Africa : " + v['gene'] :
            print("ATTENTION DIFFERENCE ====> %s:%s" % (gene, v['gene']))
            gene = "NCBI : %s | H3Africa : %s" % (gene, v['gene'])
        if ch != v['ch']:
            print("ATTENTION DIFFERENCE ====> %s:%s" % (ch, v['ch']))
            ch = "NCBI : %s | H3Africa : %s" % (ch, v['ch'])
        if pos != k:
            print("ATTENTION DIFFERENCE ====> %s:%s" % (pos, k))
            pos = "NCBI : %s | H3Africa : %s" % (pos, k)
    else:
        gene = "H3Africa : " + v['gene']
        rsID = "N/A"
        ch = v['ch']
        pos = k
        decoded_ALLELE = "N/A"
    row = [gene, rsID, ch, pos, decoded_ALLELE]
    w.writerow(row)

result_file.close()
