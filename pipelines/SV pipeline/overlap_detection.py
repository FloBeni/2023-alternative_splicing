
import pandas as pd
import sys
from collections import Counter
import re


print(sys.argv)

by_intron_cds_path = sys.argv[1]
overlap_path = sys.argv[2]
fpkm_cov_path = sys.argv[3]
minor_introns_path = sys.argv[4]

# by_intron_cds_path = "/home/fbenitiere/data/Projet-SplicedVariants/Analyses/Drosophila_melanogaster/by_intron_cds.tab"
# fpkm_cov_path = "/home/fbenitiere/data//Projet-SplicedVariants//Analyses/Macaca_mulatta/by_gene_analysis.tab"
# minor_introns_path = "/home/fbenitiere/data//Projet-SplicedVariants/Analyses/Macaca_mulatta/by_minor_intron.tab"


intron_cds_table = pd.read_table(by_intron_cds_path , dtype={'into_cds': 'str'} , comment="#")
Counter(intron_cds_table.into_cds)

intron_cds_table = intron_cds_table.loc[intron_cds_table.intron_class == "major"]
intron_cds_table["start"] = intron_cds_table[["splice5","splice3"]].min(axis=1)
intron_cds_table["end"] = intron_cds_table[["splice5","splice3"]].max(axis=1)

def detection_overlap(x):
    nb_overlap = len( intron_cds_table.loc[ (intron_cds_table.start < x.end) & (intron_cds_table.end > x.start) & (intron_cds_table.gene_id == x.gene_id) ] ) -1
    return(nb_overlap)


intron_cds_table["overlap"] = intron_cds_table.apply( detection_overlap , axis = 1 )
Counter(intron_cds_table.overlap)

intron_cds_table["id"] = intron_cds_table.seqname + ";"+intron_cds_table.gene_id + ";"+intron_cds_table.strand.astype(str)  + ";"+intron_cds_table.splice5.astype(str)  + ";"+intron_cds_table.splice3.astype(str)


fpkm_cov = pd.read_table(fpkm_cov_path , comment="#")
fpkm_cov = fpkm_cov.loc[(fpkm_cov.type == "gene" ) &  [re.search('gene_biotype=protein_coding', x) for x in fpkm_cov.attributes] ]
fpkm_cov.set_index(fpkm_cov.gene_id  , inplace=True)

candidate = intron_cds_table.loc[(intron_cds_table.intron_class == "major") & (intron_cds_table.into_cds == "True") & intron_cds_table.gene_id.isin(fpkm_cov.gene_id)]

candidate = candidate.loc[candidate.splice_variant_rate > 0.05]

minor_introns = pd.read_table(minor_introns_path , comment="#")
minor_introns = minor_introns.loc[(minor_introns.intron_class == "minor" ) & minor_introns.gene_id.isin(candidate.gene_id) ]

def major_id(intron):
    if intron.which_shared_site == "splice5":
      id_maj = intron.seqname + ";"+intron.gene_id + ";"+ str(int(intron.strand))  + ";" + str(int(intron.splice5))  + ";"+str(int(intron.coord_specific_major_site))
    elif intron.which_shared_site == "splice3":
      id_maj = intron.seqname + ";"+intron.gene_id + ";"+str(int(intron.strand))  + ";"+str(int(intron.coord_specific_major_site))  + ";"+str(int(intron.splice3))
    else :
      return("both")
    return(id_maj)


minor_introns['major_correspond'] = [major_id(intron) for intron in minor_introns.itertuples()]

minor_introns = minor_introns.loc[(minor_introns.which_shared_site == "both" ) | minor_introns.major_correspond.isin(candidate.id) ]

def abundant_sv(intron):
  minor_candidate = minor_introns.loc[(minor_introns.seqname == intron.seqname) & (minor_introns.gene_id == intron.gene_id) & ((minor_introns.splice5 == int(intron.splice5)) | (minor_introns.splice3 == int(intron.splice3)))].copy()
  minor_candidate.mira = minor_candidate.n1 / (int(intron.n1) + int(intron.n2_spl3) + int(intron.n2_spl5))
  return(sum(minor_candidate.mira > 0.05))


candidate['abundant_sv'] = [abundant_sv(intron) for intron in candidate.itertuples()]

Counter(candidate.abundant_sv == 0)

candidate = candidate.loc[candidate.abundant_sv != 0]

intron_cds_table["have_abundant_sv"] = intron_cds_table.id.isin(candidate.id)
Counter(intron_cds_table.have_abundant_sv )

intron_cds_table.to_csv(overlap_path, index=False, sep="\t")
