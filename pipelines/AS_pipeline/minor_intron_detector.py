

import pandas as pd
import sys
import statistics
from collections import Counter



print(sys.argv)

intron_cds_path = sys.argv[1]
minor_intron_path = sys.argv[2]


intron_table = pd.read_table(intron_cds_path, dtype={'into_cds': 'str'} , comment="#")
intron_table["n2"] = intron_table.n2_spl3 + intron_table.n2_spl5
intron_table["start"] = intron_table[["splice5","splice3"]].min(axis=1)
intron_table["end"] = intron_table[["splice5","splice3"]].max(axis=1)

major_intron_table = intron_table.loc[((intron_table.intron_class == "major") & (intron_table.into_cds == "True"))].copy()
minor_intron_table = intron_table.loc[intron_table.intron_class == "minor"].copy()

minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice5.astype(str)  , inplace=True)

major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice5.astype(str)  , inplace=True)

minor_intron_table["share_major_spl5_site"] = False
minor_intron_table.loc[minor_intron_table.index.isin(major_intron_table.index),"share_major_spl5_site"] = True



minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice3.astype(str)  , inplace=True)

major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice3.astype(str)  , inplace=True)
minor_intron_table["share_major_spl3_site"] = False
minor_intron_table.loc[minor_intron_table.index.isin(major_intron_table.index),"share_major_spl3_site"] = True




minor_intron_table["which_shared_site"] = None
minor_intron_table.loc[ minor_intron_table.share_major_spl5_site, "which_shared_site"] ="splice5"
minor_intron_table.loc[ minor_intron_table.share_major_spl3_site , "which_shared_site"] ="splice3"
minor_intron_table.loc[ minor_intron_table.share_major_spl3_site & minor_intron_table.share_major_spl5_site , "which_shared_site"] ="both"


minor_intron_table = minor_intron_table.loc[~minor_intron_table.which_shared_site.isnull()]


minor_intron_table["coord_specific_major_site"] = None
minor_intron_table["n1_major"] = None
minor_intron_table["n2_major"] = None

major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice5.astype(str)  , inplace=True)
minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice5.astype(str)  , inplace=True)
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5","coord_specific_major_site"] = major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice5"].index].splice3.values
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5","n1_major"] = major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice5"].index].n1.values
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5","n2_major"] = major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice5"].index].n2.values


major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice3.astype(str)  , inplace=True)
minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice3.astype(str)  , inplace=True)
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","coord_specific_major_site"] = major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3"].index].splice5.values
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","n1_major"] = major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3"].index].n1.values
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","n2_major"] = major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3"].index].n2.values





minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5","distance_from_major"] = (minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5","splice3"] - minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5","coord_specific_major_site"]).abs()
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","distance_from_major"] = (minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","splice5"] - minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","coord_specific_major_site"]).abs()

minor_intron_table["frame_shift"] = minor_intron_table.distance_from_major % 3



minor_intron_table["mira"] = minor_intron_table.n1 / (minor_intron_table.n1_major + minor_intron_table.n2_major)



def detect_mira_both(intron):
	major_candidate = major_intron_table.loc[(major_intron_table.seqname == intron.seqname) & (major_intron_table.gene_id == intron.gene_id) & ((major_intron_table.splice5 == int(intron.splice5)) | (major_intron_table.splice3 == int(intron.splice3)))].copy()
	major_candidate["mira"] = int(intron.n1) / (major_candidate.n1.values + major_candidate.n2.values)
	return(statistics.mean(major_candidate.mira.values))

minor_intron_table.loc[minor_intron_table.which_shared_site == "both","mira"] = [detect_mira_both(intron) for intron in minor_intron_table.loc[minor_intron_table.which_shared_site == "both",].itertuples()]


### DETERMINER SI BORNE DANS LINTRON MAJORITAIRE ??

minor_intron_table["specific_splsite_into_shared_major"] = None

major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice3.astype(str)  , inplace=True)
minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice3.astype(str)  , inplace=True)
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice3","specific_splsite_into_shared_major"] =	((
																															 major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3"].index].start <= minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3","start"])
																													 &
																													 (major_intron_table.loc[minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3"].index].end >= minor_intron_table.loc[ minor_intron_table.which_shared_site == "splice3","end"]))

major_intron_table.set_index(
	major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
		str) + '|' + major_intron_table.splice5.astype(str), inplace=True)
minor_intron_table.set_index(
	minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
		str) + '|' + minor_intron_table.splice5.astype(str), inplace=True)
minor_intron_table.loc[minor_intron_table.which_shared_site == "splice5", "specific_splsite_into_shared_major"] = ((
																														   major_intron_table.loc[
																															   minor_intron_table.loc[
																																   minor_intron_table.which_shared_site == "splice5"].index].start <=
																														   minor_intron_table.loc[
																															   minor_intron_table.which_shared_site == "splice5", "start"])
																												   &
																												   (
																														   major_intron_table.loc[
																															   minor_intron_table.loc[
																																   minor_intron_table.which_shared_site == "splice5"].index].end >=
																														   minor_intron_table.loc[
																															   minor_intron_table.which_shared_site == "splice5", "end"]))




major_intron_table = intron_table.loc[intron_table.intron_class == "major"].copy()

major_intron_table = major_intron_table.sort_values(ascending=True,by=["gene_id","seqname" ,"strand","start"])

major_intron_table["start_next_major"] = None
major_intron_table["gene_next_major"] = None

major_intron_table.iloc[0:(len(major_intron_table)-1),major_intron_table.columns.get_loc('start_next_major')] = major_intron_table.iloc[1:len(major_intron_table)].start
major_intron_table.iloc[0:(len(major_intron_table)-1),major_intron_table.columns.get_loc('gene_next_major')] = major_intron_table.iloc[1:len(major_intron_table)].gene_id


major_intron_table["end_previous_major"] = None
major_intron_table["gene_previous_major"] = None

major_intron_table.iloc[1:len(major_intron_table),major_intron_table.columns.get_loc('end_previous_major')] = major_intron_table.iloc[0:(len(major_intron_table)-1)].end
major_intron_table.iloc[1:len(major_intron_table),major_intron_table.columns.get_loc('gene_previous_major')] = major_intron_table.iloc[0:(len(major_intron_table)-1)].gene_id



minor_intron_table["criptic_intron"] = None

major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice3.astype(str)  , inplace=True)
minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice3.astype(str)  , inplace=True)
INDEX = minor_intron_table.which_shared_site == "splice3"
minor_intron_table.loc[INDEX,"criptic_intron"] = 	(( major_intron_table.loc[minor_intron_table.loc[ INDEX ].index].start_next_major <= minor_intron_table.loc[ INDEX ,"end"] )
													 & (major_intron_table.loc[minor_intron_table.loc[ INDEX ].index].gene_next_major == minor_intron_table.loc[ INDEX,"gene_id"])).values|(
															(major_intron_table.loc[minor_intron_table.loc[ INDEX ].index].end_previous_major >= minor_intron_table.loc[ INDEX,"start"])
															& (major_intron_table.loc[minor_intron_table.loc[ INDEX ].index].gene_previous_major == minor_intron_table.loc[ INDEX,"gene_id"])).values



major_intron_table.set_index(major_intron_table.seqname + '|' + major_intron_table.gene_id + '|' + major_intron_table.strand.astype(
	str) + '|' + major_intron_table.splice5.astype(str)  , inplace=True)
minor_intron_table.set_index(minor_intron_table.seqname + '|' + minor_intron_table.gene_id + '|' + minor_intron_table.strand.astype(
	str) + '|' + minor_intron_table.splice5.astype(str)  , inplace=True)
INDEX = minor_intron_table.which_shared_site == "splice5"
minor_intron_table.loc[INDEX,"criptic_intron"] = 	((major_intron_table.loc[minor_intron_table.loc[ INDEX].index].start_next_major <= minor_intron_table.loc[ INDEX,"end"])
													 & (major_intron_table.loc[minor_intron_table.loc[ INDEX].index].gene_next_major == minor_intron_table.loc[ INDEX,"gene_id"])).values|	(
															(major_intron_table.loc[minor_intron_table.loc[ INDEX].index].end_previous_major >= minor_intron_table.loc[ INDEX,"start"])
															& (major_intron_table.loc[minor_intron_table.loc[ INDEX].index].gene_previous_major == minor_intron_table.loc[ INDEX,"gene_id"])).values


Counter(minor_intron_table.criptic_intron)


minor_intron_table = minor_intron_table[['gene_id','seqname','strand','splice5','splice3','n1','n2_spl5',
										 'n2_spl3','n3_spl3','n3_spl5','splice_variant_rate',
										 'nonsplice_variant_rate','intron_class','into_cds','n2','start','end',
										 'which_shared_site','coord_specific_major_site','n1_major','n2_major',
										 'distance_from_major','frame_shift',
										 'mira','specific_splsite_into_shared_major','criptic_intron']]

minor_intron_table.to_csv(	minor_intron_path  , index = False, sep = "\t")
