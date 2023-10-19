
import pandas as pd
import sys
from collections import Counter
import gffpandas.gffpandas as gffpd

print(sys.argv)

path_input = sys.argv[1]
list_Acc_path = sys.argv[2]
gff_path = sys.argv[3]
by_intron_db_path = sys.argv[4]
by_gene_db_path = sys.argv[5]

list_Acc = pd.read_table(list_Acc_path)
list_Acc.set_index(list_Acc.BioProjet,inplace=True)





annotation = gffpd.read_gff3(gff_path)
print(annotation.header)
print(annotation.df)
gene_annotation = annotation.filter_feature_of_type(['gene','pseudogene'])
gene_annotation = gene_annotation.attributes_to_columns()
gene_annotation.set_index( gene_annotation.ID,inplace=True )

gene_annotation = gene_annotation.loc[~gene_annotation.ID.duplicated()]



passage = 0

for rna_samp in list_Acc.SRA_accession_ID:
	print(rna_samp)
	intron_samp_table = pd.read_table(path_input + rna_samp + "/n1_n2_n3_by_intron.tab")
	gene_samp_table = pd.read_table(path_input + rna_samp + "/coverage_fpkm_by_genes.tab")
	if passage == 0:
		by_intron_compilation = intron_samp_table
		by_genes_compilation = gene_samp_table

		by_genes_compilation["attributes"] = gene_annotation.loc[by_genes_compilation.gene_id]["attributes"].values
		by_genes_compilation["seq_id"] = gene_annotation.loc[by_genes_compilation.gene_id]["seq_id"].values
		by_genes_compilation["start"] = gene_annotation.loc[by_genes_compilation.gene_id]["start"].values
		by_genes_compilation["end"] = gene_annotation.loc[by_genes_compilation.gene_id]["end"].values
		by_genes_compilation["strand"] = gene_annotation.loc[by_genes_compilation.gene_id]["strand"].values
		by_genes_compilation["type"] = gene_annotation.loc[by_genes_compilation.gene_id]["type"].values
		by_genes_compilation = by_genes_compilation[['gene_id', 'gene_name','seq_id','start','end','strand','type','attributes','exon_coverage']]

		by_genes_compilation["exon_coverage_" + rna_samp] = gene_samp_table.exon_coverage
		by_genes_compilation["fpkm_" + rna_samp] = gene_samp_table.fpkm


		by_intron_compilation['sum_n1'] =  intron_samp_table.n1
		by_intron_compilation['sum_n2'] =  intron_samp_table.n2_spl5 + intron_samp_table.n2_spl3
		by_intron_compilation['sum_n3'] =  intron_samp_table.n3_spl5 + intron_samp_table.n3_spl3


		by_intron_compilation = by_intron_compilation[["gene_id",	"seqname","strand","splice5","splice3",	"sum_n1", "sum_n2",	"sum_n3"]]

		by_intron_compilation["n1_" + rna_samp] = intron_samp_table.n1
		by_intron_compilation["n2_" + rna_samp] = intron_samp_table.n2_spl5 + intron_samp_table.n2_spl3
		by_intron_compilation["n3_" + rna_samp] = intron_samp_table.n3_spl5 + intron_samp_table.n3_spl3

		fpkm = pd.DataFrame(gene_samp_table.fpkm.values , columns = [rna_samp] )
		passage = 1
	else:
		by_intron_compilation.sum_n1 = by_intron_compilation.sum_n1  + intron_samp_table.n1
		by_intron_compilation.sum_n2 = by_intron_compilation.sum_n2  + intron_samp_table.n2_spl5 + intron_samp_table.n2_spl3
		by_intron_compilation.sum_n3 = by_intron_compilation.sum_n3  + intron_samp_table.n3_spl5 + intron_samp_table.n3_spl3

		by_intron_compilation["n1_" + rna_samp] = intron_samp_table.n1
		by_intron_compilation["n2_" + rna_samp] = intron_samp_table.n2_spl5 + intron_samp_table.n2_spl3
		by_intron_compilation["n3_" + rna_samp] = intron_samp_table.n3_spl5 + intron_samp_table.n3_spl3

		by_genes_compilation["exon_coverage_" + rna_samp] = gene_samp_table.exon_coverage
		by_genes_compilation["fpkm_" + rna_samp] = gene_samp_table.fpkm
		by_genes_compilation.exon_coverage = by_genes_compilation.exon_coverage + gene_samp_table.exon_coverage

		fpkm[ rna_samp ] = gene_samp_table.fpkm.values

by_genes_compilation["median_fpkm"] = fpkm.median(axis=1)
by_genes_compilation["mean_fpkm"] = fpkm.mean(axis=1)

by_intron_compilation["splice_variant_rate"] = ( by_intron_compilation.sum_n2 ) / (by_intron_compilation.sum_n2 + by_intron_compilation.sum_n1)
by_intron_compilation.loc[(by_intron_compilation.sum_n2 + by_intron_compilation.sum_n1) < 10,"splice_variant_rate"] = None
by_intron_compilation["nonsplice_variant_rate"] = ( by_intron_compilation.sum_n3 ) / (by_intron_compilation.sum_n3 + 2*by_intron_compilation.sum_n1)
by_intron_compilation.loc[(by_intron_compilation.sum_n3/2 + by_intron_compilation.sum_n1) < 10,"nonsplice_variant_rate"] = None

by_intron_compilation["intron_class"] = "Unclassified"
by_intron_compilation.loc[((by_intron_compilation.sum_n1 > 0) & (by_intron_compilation.nonsplice_variant_rate < 0.5) &
						   (by_intron_compilation.splice_variant_rate < 0.5)), "intron_class"] = "major"

by_intron_compilation.loc[((by_intron_compilation.sum_n1 > 0) & (by_intron_compilation.nonsplice_variant_rate >= 0.5) |
						   (by_intron_compilation.splice_variant_rate >= 0.5)), "intron_class"] = "minor"



with open(by_intron_db_path, 'w') as f:
    f.write(annotation.header + "## RNAseq SRA list: " +  " ".join([str(elem) for elem in list_Acc.SRA_accession_ID]) +"\n"+
			"## n1: reads containing the spliced intron\n"+
			"## n2: reads containing variants\n"+
			"## n3: reads containing retention +10 / -10 bp\n"+
			"## splice_variant_rate: n2 / (n1 +n2)\n"+
			"## nonsplice_variant_rate: n3 / (n1*2 +n3)\n"+
			"## intron_class: major n1 > n2 and 2*n1 > n3\n"
			)
    f.close

by_intron_compilation.to_csv(by_intron_db_path , mode="a", index=False, sep="\t")


with open(by_gene_db_path, 'w') as f:
    f.write(annotation.header + "## RNAseq SRA list: " +  " ".join([str(elem) for elem in list_Acc.SRA_accession_ID]) +"\n"+
			"## median_fpkm: median FPKM across RNAseq\n"+
			"## mean_fpkm: mean FPKM across RNAseq\n"+
			"## std_fpkm: standard deviation FPKM across RNAseq\n"+
			"## exon_coverage: exonic read per bp (all RNAseq pooled)\n"
			)
    f.close

by_genes_compilation.to_csv(by_gene_db_path,mode="a", index=False, sep="\t")
