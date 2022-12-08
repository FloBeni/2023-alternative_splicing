import gffpandas.gffpandas as gffpd
import pandas as pd
import sys
from collections import Counter
import Bio.SeqIO as IO
from itertools import product
from itertools import compress
from os.path import isfile, join
import gffpandas.gffpandas as gffpd
import numpy as np

print( sys.argv )
genome_path = sys.argv[1]
intronCoord_path = sys.argv[2]
gff_path = sys.argv[3]
trinucl_intronic_path = sys.argv[4]


annotation = gffpd.read_gff3(gff_path)
print(annotation.header)
print(annotation.df)

gene_annotation = annotation.filter_feature_of_type(['CDS'])
gene_annotation = gene_annotation.attributes_to_columns()
gene_annotation.set_index( gene_annotation.Parent , inplace=True )

gene_annotation = gene_annotation.loc[~gene_annotation.protein_id.duplicated()]


# gene_annotation = gene_annotation.loc[np.invert(gene_annotation.Parent.values == None)] # Sil ny a pas de PARENT on ne garde pas la ligne...peut etre lie hastag
gene_annotation = gene_annotation.loc[ [ ";Parent=" in i for i in gene_annotation.attributes]] # Sil ny a pas de PARENT on ne garde pas la ligne...peut etre lie hastag
gene_annotation = gene_annotation.loc[ [ ";protein_id=" in i for i in gene_annotation.attributes] ]

if any(gene_annotation.protein_id.values == None):
    print(Counter(gene_annotation.protein_id.values == None))
    print(gene_annotation.loc[gene_annotation.protein_id.values == None])
    quit()

intron_library = pd.read_table(intronCoord_path)

intron_id = [ i + ':CDS' for i in gene_annotation.Parent]
intron_identified_to_prot = [ any([ j in i for j in intron_id ]) for i in intron_library.Transcripts ]

intron_library = intron_library.loc[ intron_identified_to_prot ]


def detect_protein(x):
    list_prot = list(compress(intron_id, [ i in x.Transcripts for i in intron_id]))
    return( ";".join(gene_annotation.loc[[i.replace(':CDS', '') for i in list_prot ]].protein_id) )

intron_library["protein_id"] = intron_library.apply(detect_protein,axis=1)

intron_library["splice5"] = intron_library.Start
intron_library["splice3"] = intron_library.End
intron_library.loc[intron_library.Strand == -1,"splice5"] = intron_library.loc[intron_library.Strand == -1].End
intron_library.loc[intron_library.Strand == -1,"splice3"] = intron_library.loc[intron_library.Strand == -1].Start

intron_library = intron_library[['Genes','Chr','Strand','splice5','splice3','Transcripts','protein_id']]

intron_library.rename(columns={'Genes': 'gene_id', 'Chr': 'seqname','Strand': 'strand','Transcripts': 'transcript_id'}, inplace=True)


li = ['A', 'T', 'G' , 'C']
combi = ["combi"]
for comb in product(li, repeat=3):
    combi.append(''.join(comb))
combi = combi[1:66]


genome = IO.to_dict(IO.parse(genome_path, "fasta"))

triplet_freq = pd.DataFrame()

for intron in intron_library.itertuples():
    print(intron.strand)
    if intron.strand == 1:
        intron_seq = genome[intron.seqname].seq[(intron.splice5+1):(intron.splice3-2)].upper()
        print(intron_seq)
    else:
        intron_seq = str(genome[intron.seqname].seq[(intron.splice3+1):(intron.splice5-2)][::-1].upper())
        intron_seq = intron_seq.translate(intron_seq.maketrans("ATGC","TACG"))
        print(intron_seq)
    occurences = Counter([intron_seq[i:i + 3].upper() for i in range(0, len(intron_seq)-2, 1)])
    nucl_occu = Counter([intron_seq[i].upper() for i in range(0, len(intron_seq), 1)])
    freq_triplet = [occurences[triplet] for triplet in combi]
    freq_nucl = [nucl_occu[nucl] for nucl in li]
    data_prot = pd.DataFrame([ list(intron[1:]) + freq_nucl + freq_triplet ],
                             columns = ['gene_id', 'seqname', 'strand', 'splice5', 'splice3', 'transcript_id','protein_id','Ai', 'Ti', 'Gi', 'Ci'] + combi )
    triplet_freq = triplet_freq.append(data_prot)

triplet_freq.to_csv(trinucl_intronic_path, index=False, sep="\t")
