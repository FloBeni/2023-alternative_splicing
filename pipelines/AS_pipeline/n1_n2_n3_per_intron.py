import pandas as pd
import sys
from collections import Counter

print(sys.argv)
nsvr_path = sys.argv[1]
svr_path = sys.argv[2]
output_path = sys.argv[3]

# nsvr_path = "/home/florianbenitiere/data/Projet-SplicedVariants/Analyses-RNAseq/Homo_sapiens/ERR2598059/NSVR_Frequency_inclusive.txt"
# svr_path = "/home/florianbenitiere/data/Projet-SplicedVariants/Analyses-RNAseq/Homo_sapiens/ERR2598059/SVR_Frequency_inclusive.txt"

svr_table = pd.read_table(svr_path)
svr_table.columns = ["gene_id","seqname","strand","splice5","splice3","n1","otherspl5coord","otherspl3coord","n2_spl5","n2_spl3"]
svr_table = svr_table.loc[svr_table.gene_id.astype(str) != "nan"]

svr_table.set_index(svr_table.seqname + '|' + svr_table.gene_id + '|' + svr_table.splice5.astype(str) + '|' + svr_table.splice3.astype(str)  , inplace=True) # pas le strand car gene15842 ches droso bug -1 et intron +1


nsvr_table = pd.read_table(nsvr_path)
nsvr_table.columns = ["gene_id","seqname","start","end","strand","n3_spl5","n3_spl3"]
nsvr_table = nsvr_table.loc[nsvr_table.gene_id.astype(str) != "nan"]
nsvr_table = nsvr_table.assign(gene_id = nsvr_table['gene_id'].str.split(',')).explode( 'gene_id' )
nsvr_table = nsvr_table.loc[nsvr_table.gene_id.astype(str) != "NA"]

nsvr_table["splice5"] = None
nsvr_table["splice3"] = None
nsvr_table.loc[nsvr_table.strand == 1,"splice5"] = nsvr_table.loc[nsvr_table.strand == 1,"start"].values
nsvr_table.loc[nsvr_table.strand == -1,"splice5"] = nsvr_table.loc[nsvr_table.strand == -1,"end"].values
nsvr_table.loc[nsvr_table.strand == 1,"splice3"] = nsvr_table.loc[nsvr_table.strand == 1,"end"].values
nsvr_table.loc[nsvr_table.strand == -1,"splice3"] = nsvr_table.loc[nsvr_table.strand == -1,"start"].values

nsvr_table.set_index(nsvr_table.seqname + '|' + nsvr_table.gene_id + '|' + nsvr_table.splice5.astype(str) + '|' + nsvr_table.splice3.astype(str)  , inplace=True)

svr_table[ "n3_spl3" ] = nsvr_table.loc[svr_table.index, "n3_spl3" ]
svr_table[ "n3_spl5" ] = nsvr_table.loc[svr_table.index, "n3_spl5" ]


svr_table = svr_table[[ 'gene_id' , 'seqname' , 'strand' , 'splice5', 'splice3' , 'n1' , 'n2_spl5' , 'n2_spl3' , 'n3_spl5' , 'n3_spl3' ]]
svr_table.to_csv(output_path, index=False, sep="\t")
