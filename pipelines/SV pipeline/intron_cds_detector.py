import gffpandas.gffpandas as gffpd
import pandas as pd
import sys
from collections import Counter

print(sys.argv)
gff_path = sys.argv[1]
by_intron_analysis_path = sys.argv[2]
list_Acc_path = sys.argv[3]
gene_coord_path = sys.argv[4]
by_intron_cds_detection_path = sys.argv[5]

annotation = gffpd.read_gff3(gff_path)
print(annotation.header)
print(annotation.df)

cds_annotation = annotation.filter_feature_of_type(['CDS'])
cds_annotation = cds_annotation.attributes_to_columns()

mRNA_annotation = annotation.filter_feature_of_type(['mRNA'])
mRNA_annotation = mRNA_annotation.attributes_to_columns()

cds_annotation["corresponding_gene"] = cds_annotation.Parent
if len(mRNA_annotation) != 0:
    mRNA_annotation.set_index( mRNA_annotation.ID , inplace=True )
    cds_annotation.loc[ cds_annotation.Parent.isin(mRNA_annotation.index),"corresponding_gene" ] = mRNA_annotation.loc[ cds_annotation.loc[ cds_annotation.Parent.isin(mRNA_annotation.index),"corresponding_gene"] ,"Parent" ].values


cds_annotation.set_index(cds_annotation.corresponding_gene)

gene_df = cds_annotation.groupby(['corresponding_gene']).start.min().to_frame()
gene_df["end"] = cds_annotation.groupby(['corresponding_gene']).end.max()

gene_info = cds_annotation.loc[ ~cds_annotation.corresponding_gene.duplicated() ]
gene_info.set_index( gene_info.corresponding_gene,inplace=True )

gene_df["strand"] = gene_info.loc[gene_df.index,"strand"]
gene_df["seq_id"] = gene_info.loc[gene_df.index,"seq_id"]


gene_df["mini_start_codon"] = 0
gene_df["maxi_stop_codon"] = 0
gene_df.loc[gene_df.strand == "+","mini_start_codon"] = gene_df.loc[gene_df.strand == "+","start"].values
gene_df.loc[gene_df.strand == "-","mini_start_codon"] = gene_df.loc[gene_df.strand == "-","end"].values
gene_df.loc[gene_df.strand == "+","maxi_stop_codon"] = gene_df.loc[gene_df.strand == "+","end"].values
gene_df.loc[gene_df.strand == "-","maxi_stop_codon"] = gene_df.loc[gene_df.strand == "-","start"].values
Counter(gene_df.mini_start_codon == 0)
Counter(gene_df.maxi_stop_codon == 0)


gene_df.to_csv(gene_coord_path, index=True, sep="\t")

intron_metazoa = pd.read_table(by_intron_analysis_path)
Gene_with_CDS = intron_metazoa.gene_id.isin(gene_df.index)

intron_metazoa["gene_start"] = None
intron_metazoa.loc[Gene_with_CDS , "gene_start" ] = gene_df.loc[ intron_metazoa.loc[Gene_with_CDS , "gene_id" ] , "start" ].values
intron_metazoa.loc[Gene_with_CDS , "gene_end" ] = gene_df.loc[ intron_metazoa.loc[Gene_with_CDS , "gene_id" ] , "end" ].values

intron_metazoa["start"] = intron_metazoa[["splice5","splice3"]].min(axis=1)
intron_metazoa["end"] = intron_metazoa[["splice5","splice3"]].max(axis=1)

intron_metazoa["sup_start"] = None
intron_metazoa.loc[Gene_with_CDS,"sup_start"] = intron_metazoa.loc[Gene_with_CDS].start >= intron_metazoa.loc[Gene_with_CDS].gene_start
intron_metazoa.loc[Gene_with_CDS,"inf_end"] = intron_metazoa.loc[Gene_with_CDS].end <= intron_metazoa.loc[Gene_with_CDS].gene_end

intron_metazoa["into_cds"] = "undetermined"
intron_metazoa.loc[Gene_with_CDS,"into_cds"] = intron_metazoa.loc[Gene_with_CDS,"inf_end"] & intron_metazoa.loc[Gene_with_CDS,"sup_start"]

Counter(intron_metazoa.into_cds)
intron_metazoa = intron_metazoa[["gene_id",	"seqname","strand","splice5","splice3",	"n1",
                                 "n2_spl5",	"n2_spl3","n3_spl3","n3_spl5","splice_variant_rate",
                                 "nonsplice_variant_rate","intron_class","into_cds"]]




list_Acc = pd.read_table(list_Acc_path)
list_Acc.set_index(list_Acc.BioProjet,inplace=True)

with open(by_intron_cds_detection_path, 'w') as f:
    f.write(annotation.header + "## RNAseq SRA list: " +  " ".join([str(elem) for elem in list_Acc.SRA_accession_ID]) +"\n"+
			"## n1: reads containing the spliced intron\n"+
			"## n2: reads containing variants\n"+
			"## n3: reads containing retention +10 / -10 bp\n"+
			"## splice_variant_rate: n2 / (n1 +n2)\n"+
			"## nonsplice_variant_rate: n3 / (n1*2 +n3)\n"+
			"## intron_class: major n1 > n2 and 2*n1 > n3\n"
			)
    f.close

intron_metazoa.loc[intron_metazoa.nonsplice_variant_rate.isna(), 'nonsplice_variant_rate'] = "NA"
intron_metazoa.loc[intron_metazoa.splice_variant_rate.isna(), 'splice_variant_rate'] = "NA"

intron_metazoa.to_csv(by_intron_cds_detection_path,mode='a', index=False, sep="\t")
