import gffpandas.gffpandas as gffpd
import pandas as pd
import sys
from collections import Counter
import Bio.SeqIO as IO
from itertools import product
from os.path import isfile, join

print(sys.argv)
gff_path = sys.argv[1]
tRNAscan_path = sys.argv[2]

annotation = gffpd.read_gff3(gff_path)
print(annotation.header)

tRNA_annotation = annotation.filter_feature_of_type(['tRNA'])
tRNA_annotation = tRNA_annotation.attributes_to_columns()
print(tRNA_annotation.head)
tRNA_annotation = tRNA_annotation.loc[tRNA_annotation.source == "tRNAscan-SE"]

if len(tRNA_annotation) == 0:
    print(tRNAscan_path)
    open(tRNAscan_path, 'a').close()
else:
    tRNA_annotation["gene_id"] = tRNA_annotation.Parent
    tRNA_annotation = tRNA_annotation[["gene_id","start","end","seq_id","source","type","product","gene","Note"]]
    tRNA_annotation["anticodon"] = ""
    tRNA_annotation.loc[[ 'anticodon ' in tRNA for tRNA in tRNA_annotation.Note],"anticodon"] = [f.split('anticodon ')[1].split(')')[0] for f in tRNA_annotation.Note if 'anticodon ' in f]
    tRNA_annotation.loc[[ 'anticodon ' not in tRNA for tRNA in tRNA_annotation.Note],"anticodon"] = [f.split('(')[1].split(')')[0] for f in tRNA_annotation.Note if 'anticodon ' not in f]
    tRNA_annotation.to_csv(tRNAscan_path, index=False, sep="\t")


