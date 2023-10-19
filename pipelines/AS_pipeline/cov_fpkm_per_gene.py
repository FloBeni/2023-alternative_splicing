import pandas as pd
import sys
from collections import Counter

print(sys.argv)
gtf_path = sys.argv[1]
gene_intron_coverage_path = sys.argv[2]
gene_exon_coverage_path = sys.argv[3]
gene_fpkm_path = sys.argv[4]
gene_analysis_report_path = sys.argv[5]

annotation_gtf = pd.read_table(gtf_path, names=['seqname', 'source', 'feature', 'start', 'end', 'scrore', 'strand', 'frame','attribute'])

annotation_gtf = annotation_gtf.loc[ ~annotation_gtf.attribute.duplicated() ]
annotation_gtf = annotation_gtf.loc[["gene_id" in attr for attr in annotation_gtf.attribute]]
df = annotation_gtf.attribute.str.split(";", expand=True)
annotation_gtf = annotation_gtf.join(df)
# annotation_gtf = annotation_gtf.loc[~(annotation_gtf[1] == '')]
# annotation_gtf = annotation_gtf.loc[~(annotation_gtf[2] == '')]
annotation_gtf['gene_id'] = [x.replace(' gene_id "', '').replace('"', '') for x in annotation_gtf[1].values]
annotation_gtf['gene_name'] = [x.replace(' gene_name "', '').replace('"', '') for x in annotation_gtf[2].values]
annotation_gtf = annotation_gtf.loc[ ~annotation_gtf.gene_id.duplicated() ]

gene_fpkm = pd.read_table(gene_fpkm_path)
gene_fpkm = gene_fpkm.loc[~gene_fpkm.tracking_id.duplicated(keep=False)] # Supprime les duplicats dans FPKM car on ne sait pas lequel choisir
gene_fpkm.set_index(gene_fpkm.tracking_id,inplace=True)

ISIN_fpkm = annotation_gtf.gene_id.isin(gene_fpkm.tracking_id)


annotation_gtf.loc[ISIN_fpkm,"fpkm"] = gene_fpkm.loc[ annotation_gtf.loc[ISIN_fpkm,"gene_id"], "FPKM"].values

annotation_gtf.loc[ISIN_fpkm,"fpkm_conf_lo"] =  gene_fpkm.loc[ annotation_gtf.loc[ISIN_fpkm,"gene_id"], "FPKM_conf_lo"].values

annotation_gtf.loc[ISIN_fpkm,"fpkm_conf_hi"] =   gene_fpkm.loc[ annotation_gtf.loc[ISIN_fpkm,"gene_id"], "FPKM_conf_hi"].values

gene_exon_coverage = pd.read_table(gene_exon_coverage_path)
gene_exon_coverage.set_index(gene_exon_coverage.GeneID,inplace=True)
ISIN_exoncov = annotation_gtf.gene_id.isin(gene_exon_coverage.GeneID)

annotation_gtf.loc[ISIN_exoncov,"exon_coverage"] =  gene_exon_coverage.loc[ annotation_gtf.loc[ISIN_exoncov,"gene_id"], "Coverage"].values



gene_intron_coverage = pd.read_table(gene_intron_coverage_path)
gene_intron_coverage.set_index(gene_intron_coverage.GeneID,inplace=True)
ISIN_introncov = annotation_gtf.gene_id.isin(gene_intron_coverage.GeneID)

annotation_gtf.loc[ISIN_introncov,"intron_coverage"] = gene_intron_coverage.loc[ annotation_gtf.loc[ISIN_introncov,"gene_id"], "Coverage"].values

annotation_gtf = annotation_gtf[['gene_id', 'gene_name', 'fpkm', 'fpkm_conf_lo', 'fpkm_conf_hi', 'exon_coverage', 'intron_coverage']]


annotation_gtf.to_csv( gene_analysis_report_path , index = False, sep = "\t")
