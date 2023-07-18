import gffpandas.gffpandas as gffpd
import pandas as pd
import sys
from collections import Counter
import Bio.SeqIO as IO
from itertools import product
from os.path import isfile, join

print(sys.argv)
genome_path = sys.argv[1]
library_path = sys.argv[2]
trinucl_intronic_path = sys.argv[3]

li = ['A', 'T', 'G' , 'C']
combi = ["combi"]
for comb in product(li, repeat=3):
    combi.append(''.join(comb))
combi = combi[1:66]


intron_library = pd.read_table(library_path , dtype={'into_cds': 'str'} , comment="#")

genome = IO.to_dict(IO.parse(genome_path, "fasta"))

triplet_freq = pd.DataFrame()

for intron in intron_library.itertuples():
    if intron.strand == 1:
        intron_seq = genome[intron.seqname].seq[(intron.splice5+1):(intron.splice3-2)].upper()
    else:
        intron_seq = str(genome[intron.seqname].seq[(intron.splice3+1):(intron.splice5-2)][::-1].upper())
        intron_seq = intron_seq.translate(intron_seq.maketrans("ATGC","TACG"))
    occurences = Counter([intron_seq[i:i + 3].upper() for i in range(0, len(intron_seq)-2, 1)])
    nucl_occu = Counter([intron_seq[i].upper() for i in range(0, len(intron_seq), 1)])
    freq_triplet = [occurences[triplet] for triplet in combi]
    freq_nucl = [nucl_occu[nucl] for nucl in li]
    data_prot = pd.DataFrame([ list(intron[1:]) + freq_nucl + freq_triplet ],
                             columns = ['gene_id', 'seqname', 'strand', 'splice5', 'splice3', 'n1', 'n2_spl5',
       'n2_spl3', 'n3_spl3', 'n3_spl5', 'splice_variant_rate',
       'nonsplice_variant_rate', 'intron_class', 'into_cds', 'start', 'end',
       'overlap','Ai', 'Ti', 'Gi', 'Ci'] + combi )
    triplet_freq = triplet_freq.append(data_prot)

triplet_freq.to_csv(trinucl_intronic_path, index=False, sep="\t")
