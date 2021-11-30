import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from statistics import mean
from scipy.stats import skewnorm, norm
import numpy as np
import pandas as pd
from scipy import spatial
from pymol_test import freq_ptm_index_gen_batch_v2, fasta_reader,peptide_counting
from glob import glob
from commons import get_unique_peptide


def skew(x,e=0,w=1,a=0):
    t = (x-e) / w
    return 2 * norm.pdf(t) * norm.cdf(a*t)


# plddt_one_d = pickle.load(open('D:/data/alphafold_pdb/pLDDT_human_1d.pkl','rb'))
# sns.kdeplot(plddt_one_d,linewidth=2, alpha=.5)
# plt.show()

### plot pLDDT distribution
"""
fig, ax = plt.subplots(1, 1)


plddt_2d = pickle.load(open('D:/data/alphafold_pdb/pLDDT_human_2d.pkl','rb'))
plddt_ave = [mean(each) for each in plddt_2d]
sns.kdeplot(plddt_ave,linewidth=5, alpha=.5, color='r',ax=ax, label='average pLDDT from Alphafold2-human')

e = 95 # location
w = 55 # scale
x = np.linspace(0,100,1000)


p = skew(x,e,w,5)
ax.plot(x,p, lw=5,alpha=.5,color='k', label='optimal human average pLDDT distribution')
ax.legend()
# plt.xlim([0, 100])
plt.show()
"""

### compare distance matrix between replicates
"""

fasta_dict = fasta_reader('D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta')

df_RN = pd.read_excel('D:/data/native_protein_digestion/11182021/search_result_RN/cov_distance_each_unique_RN.xlsx',index_col=0)
df_XS = pd.read_excel('D:/data/native_protein_digestion/11182021/search_result_XS/cov_distance_each_unique_XS.xlsx',index_col=0)

protein_list_RN,protein_list_XS = df_RN.index.tolist(), df_XS.index.tolist()
print (len(protein_list_XS),len(protein_list_RN))
overlaped_protein_list = [each for each in protein_list_XS if each in protein_list_RN]
protein_dict = {each:fasta_dict[each] for each in overlaped_protein_list}

# fill null value with previous value
# df_RN_filled = df_RN.T.ffill().bfill().T
# df_XS_filled = df_XS.T.ffill().bfill().T

# forward filling and fill the rest with 100
df_RN_filled = df_RN.T.ffill().fillna(100).T
df_XS_filled = df_XS.T.ffill().fillna(100).T


peptide_tsv_RN, peptide_tsv_XS = glob('D:/data/native_protein_digestion/11182021/search_result_RN/*/peptide.tsv'), \
                                 glob('D:/data/native_protein_digestion/11182021/search_result_XS/*/peptide.tsv')
columns = ['freq_array_cossim_'+each.split('\\')[-2].split('_')[0] for each in peptide_tsv_RN]+['distance_cossim', 'Euclidean']
df_new = pd.DataFrame(index=overlaped_protein_list,columns=columns)


unique_peptide_RN,unique_peptide_XS = get_unique_peptide(peptide_tsv_RN),get_unique_peptide(peptide_tsv_XS)

for rn, xs in zip(peptide_tsv_RN,peptide_tsv_XS):
    file_name = rn.split('\\')[-2].split('_')[0]
    print (file_name)
    peptide_rn, peptide_xs = unique_peptide_RN[rn.split('\\')[-2]],unique_peptide_XS[xs.split('\\')[-2]]
    freq_array_dict_rn, freq_array_dict_xs = freq_ptm_index_gen_batch_v2(peptide_rn,protein_dict)[0], \
                                             freq_ptm_index_gen_batch_v2(peptide_xs,protein_dict)[0]
    for each in protein_dict:

        df_new.at[each,'freq_array_cossim_'+file_name] = 1-spatial.distance.cosine(freq_array_dict_rn[each],freq_array_dict_xs[each])

### compare distance cosine similarity
for each in protein_dict:
    dist_RN = df_RN_filled.loc[each,:].to_numpy()
    dist_XS = df_XS_filled.loc[each,:].to_numpy()
    df_new.at[each,'distance_cossim'] = 1-spatial.distance.cosine(dist_RN,dist_XS)
    df_new.at[each,'Euclidean'] = np.linalg.norm(dist_RN-dist_XS)
df_new.to_excel('D:/data/native_protein_digestion/11182021/cos_sim_xs_rn_2.xlsx')
"""
### add plddt column
# df = pd.read_excel('D:/data/native_protein_digestion/11182021/cos_sim_xs_rn_2.xlsx',index_col=0)
# protein_list = df.index.tolist()
# import pickle
# plddt_dict = pickle.load(open('D:/data/alphafold_pdb/pLDDT_human_dict.pkl','rb'))
#
# for prot in protein_list:
#     if prot in plddt_dict:
#         df.loc[prot,'ave_plddt'] = np.mean(plddt_dict[prot])
# df.to_excel('D:/data/native_protein_digestion/11182021/cos_sim_plddt_xs_rn.xlsx')

from pdb_operation import residue_density_cal
pdb_file = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-P61604-F1-model_v1.pdb'
residue_density_cal(pdb_file)

