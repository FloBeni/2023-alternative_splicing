import gffpandas.gffpandas as gffpd
import pandas as pd
import sys
from collections import Counter

print(sys.argv)
gff_path = sys.argv[1]
busco_csv_path = sys.argv[2]
cds_from_genomic_path = sys.argv[3]
busco_to_gene_path = sys.argv[4]

annotation = gffpd.read_gff3(gff_path)
print(annotation.header)
print(annotation.df)

cds_annotation = annotation.filter_feature_of_type(['CDS'])
cds_annotation = cds_annotation.attributes_to_columns()

mRNA_annotation = annotation.filter_feature_of_type(['mRNA'])
mRNA_annotation = mRNA_annotation.attributes_to_columns()

cds_annotation["corresponding_gene"] = cds_annotation.Parent
if len(mRNA_annotation) != 0:
    mRNA_annotation.set_index( mRNA_annotation.ID,inplace=True )
    cds_annotation.loc[ cds_annotation.Parent.isin(mRNA_annotation.index),"corresponding_gene" ] = mRNA_annotation.loc[ cds_annotation.loc[ cds_annotation.Parent.isin(mRNA_annotation.index),"corresponding_gene"] ,"Parent" ].values
cds_annotation = cds_annotation.loc[~cds_annotation.protein_id.duplicated()]
cds_annotation.set_index( cds_annotation.protein_id,inplace=True )

busco_csv = pd.read_table(busco_csv_path)
busco_csv = pd.DataFrame(busco_csv.iloc[:,0].str.split(',').tolist())
busco_csv.columns = ['busco_id', 'status', 'sequence', 'score', 'length']
busco_csv = busco_csv.iloc[4:,:]

busco_csv = busco_csv.loc[busco_csv.status != "Missing",]

busco_csv["gene_id"] = cds_annotation.loc[busco_csv.sequence , "corresponding_gene"].values

busco_csv["duplicated_id"] = busco_csv.busco_id + "|" + busco_csv.gene_id

Counter(busco_csv.duplicated_id.duplicated())

len(busco_csv.busco_id.unique())
import Bio.SeqIO as IO
record_dict = IO.to_dict(IO.parse(cds_from_genomic_path, "fasta"))
protein_id_list = []
length_cds_list = []
for key in record_dict.items():
    if "protein_id=" in key[1].description:
        protein_id = key[1].description.split("protein_id=")[1].split("]")[0]
        length_cds = len(key[1].seq)
        length_cds_list.append(length_cds)
        protein_id_list.append(protein_id)

cds_table = pd.DataFrame(length_cds_list,protein_id_list)
cds_table.columns = ["length_cds_list"]

cds_table = cds_table.sort_values(ascending=False,by=["length_cds_list"])
cds_table = cds_table.loc[~cds_table.index.duplicated()] # elimine les proteines ayant plus d'une cds en gardant la plus grande


busco_csv["cds_length"] = cds_table.loc[busco_csv.sequence,].values


busco_csv = busco_csv.sort_values(ascending=False,by=['cds_length'])
busco_csv = busco_csv.loc[~busco_csv.duplicated_id.duplicated()]


print('Counter(busco_csv.busco_id.duplicated()')
print(Counter(busco_csv.busco_id.duplicated()))
print('Counter(busco_csv.gene_id.duplicated()')
print(Counter(busco_csv.gene_id.duplicated()))


# busco_csv = busco_csv.loc[~busco_csv.busco_id.duplicated(keep=False) & ~busco_csv.gene_id.duplicated(keep=False)]

busco_csv.to_csv(busco_to_gene_path, index=False, sep="\t")
