"""
some analysis from PDB, e.g. Alphafold PLDDT,
"""

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


def map_k_r(psm_list, protein_dict):
    """
    map the start and end of each peptide
    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """

    from collections import Counter
    from collections import defaultdict
    from heapq import heappush, heappop
    import time
    import commons
    from pymol_test import automaton_matching, automaton_trie

    id_kr_mapp_dict = {}

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    separtor_pos_array = commons.separator_pos(seq_line)
    print('aho mapping...')
    aho_result = automaton_matching(automaton_trie(psm_list), seq_line)
    print('aho mapping done.')
    for tp in aho_result:
        # matched_pep = tp[2]  # without ptm site
        zero_line[tp[0] - 1] += 1
        zero_line[tp[1]] += 1

    time_start = time.time()
    for i in range(len(separtor_pos_array) - 1):
        zero_line_slice = zero_line[separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
        if np.count_nonzero(zero_line_slice) != 0:
            id_kr_mapp_dict[id_list[i]] = zero_line_slice

    return id_kr_mapp_dict


def kr_calculate(id_kr_mapp_dict, protein_dict):
    from collections import Counter
    id_kr_count_dict = {}
    k_sum, r_sum = 0, 0
    for prot in id_kr_mapp_dict:
        kr_index = np.nonzero(id_kr_mapp_dict[prot])[0]
        prot_seq = protein_dict[prot]
        kr_count_dict = Counter([prot_seq[idx] for idx in kr_index])
        k_sum += kr_count_dict['K']
        r_sum += kr_count_dict['R']
        id_kr_count_dict[prot] = {'K': kr_count_dict['K'], 'R': kr_count_dict['R']}
    return id_kr_count_dict, k_sum, r_sum


def skew(x,e=0,w=1,a=0):
    t = (x-e) / w
    return 2 * norm.pdf(t) * norm.cdf(a*t)


if __name__ == '__main__':
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

    # from pdb_operation import residue_density_cal
    # pdb_file = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-P61604-F1-model_v1.pdb'
    # residue_density_cal(pdb_file)

    ### download pdb files from pdb archive
    """
    import wget
    f = open('D:/data/pdb/human_pdb_ids.txt','r')
    pdb_id_list = f.read().split(',')
    pdb_list_len = len(pdb_id_list)
    f.close()
    
    count = 0
    for pdb in pdb_id_list:
        url = 'https://files.rcsb.org/download/'+pdb+'.pdb'
        try:
            wget.download(url, out='D:/data/pdb/pdb_human_file/')
            print (f'{count} finished, {pdb_list_len-count} to go...')
        except:
            print ('error')
            continue
        count+=1
    """

    ### read sequence from pdb fasta into dict
    """
    from pdb_operation import read_pdb_fasta
    from glob import glob
    import pickle
    pdb_seq_dict = {}
    pdb_fastas = glob('D:/data/pdb/pdb_human_sequence_file/*.fasta')
    for each in pdb_fastas:
        pdb_id = each.split('\\')[-1].split('.fasta')[0].split('_')[-1]
        seq = read_pdb_fasta(each)
        pdb_seq_dict[pdb_id] = seq
    pickle.dump(pdb_seq_dict,open('D:/data/pdb/human_pdb_seq_dict.p','wb'))
    """

    ### map pdb entries with alphafold entries
    """
    from pymol_test import pdb_file_reader
    alpha_fold_pdbs = glob('D:/data/alphafold_pdb/UP000005640_9606_HUMAN/*.pdb')
    file_len = len(alpha_fold_pdbs)
    
    alpha_pdb_dict = pickle.load(open('D:/data/alphafold_pdb/human_alpha_seq_dict.p','rb'))
    
    pdb_seq_dict = pickle.load(open('D:/data/pdb/human_pdb_seq_dict.p','rb'))
    seq_pdb_dict = {pdb_seq_dict[each]:each for each in pdb_seq_dict}
    seq_alpha_dict = {alpha_pdb_dict[each]:each for each in alpha_pdb_dict}
    
    f_out_put = open('D:/data/pdb/pdb_alpha_exactly_same.txt','w',newline='\n')
    for seq in seq_pdb_dict:
        if seq in seq_alpha_dict:
            print (seq_pdb_dict[seq],seq_alpha_dict[seq])
            f_out_put.write(seq_pdb_dict[seq]+'\t'+seq_alpha_dict[seq]+'\n')
    """

    ### K R density analysis
    """
    from statannot import add_stat_annotation
    from scipy.stats import mannwhitneyu
    kr_density_dict = pickle.load(open('D:/data/alphafold_pdb/human_file_KR_density_dict.pkl','rb'))
    k_density_array = [v for each in kr_density_dict for v in kr_density_dict[each][0].values()]
    r_density_array = [v for each in kr_density_dict for v in kr_density_dict[each][1].values()]
    u,p = mannwhitneyu(k_density_array,r_density_array)
    print (u,p)
    df_plot = pd.DataFrame(dict(residue=['K']*len(k_density_array)+['R']*len(r_density_array),
                                atom_number_within_range=k_density_array+r_density_array))
    # sns.kdeplot(
    #    data=df_plot, x="atom_number_within_range", hue="residue",
    #    fill=True, common_norm=False, palette="viridis",
    #    alpha=.5, linewidth=0,
    # )
    fig, ax = plt.subplots(1,1)
    sns.violinplot(data=df_plot,x='residue',y='atom_number_within_range',palette="viridis")
    add_stat_annotation(ax,data=df_plot,x='residue',y='atom_number_within_range',
                        box_pairs=[('K','R')],test='Mann-Whitney',text_format='star',loc='outside',verbose=2, comparisons_correction=None)
    plt.show()
    """

    # peptsv_path = glob('D:/data/deep_proteome/20200915_tryp_37C_*min/peptide.tsv')

    # unique_pep_dict = get_unique_peptide(peptsv_path)
    # print ({each:len(unique_pep_dict[each]) for each in unique_pep_dict})
    # protein_dict = fasta_reader('D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta')

    ## check unique KR cleavage frequency
    # for each in peptsv_path:
    #     file_name = each.split('\\')[-2]
    #     id_kr_mapp_dict = map_k_r(unique_pep_dict[file_name],protein_dict)
    #     print (file_name)
    #     print(kr_calculate(id_kr_mapp_dict,protein_dict)[-2:])

    ### check inslico-digested peptides
    # inslico_human_peptides_dict = pickle.load(open('D:/data/native_protein_digestion/inslico_digest_human_fasta.p','rb'))
    # inslico_peptide_list = [pep for v in inslico_human_peptides_dict.values() for pep in v]
    # print (len(inslico_peptide_list))
    # id_kr_mapp_dict = map_k_r(inslico_peptide_list,protein_dict)
    # print(kr_calculate(id_kr_mapp_dict,protein_dict)[-2:])

    ### missed cleavage analysis
    # from commons import miss_cleavage_identify
    # for each in peptsv_path:
    #     file_name = each.split('\\')[-2]
    #     kmiss_sum,r_miss_sum = miss_cleavage_identify(unique_pep_dict[file_name],regex_pattern={'K':r'K(?=[^P])','R':r'R(?=[^P])'})
    #     print (file_name,kmiss_sum,r_miss_sum)

    plot_excel = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/plot.xlsx'
    df_plot = pd.read_excel(plot_excel)
    fig, ax = plt.subplots()
    ax.set_xticks([1, 2, 3, 4])
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xticklabels(['1 hour', '2 hour', '4 hour', '20 hour'])

    sns.lineplot(ax=ax, data=df_plot, x='time', y='K/R distance to centroid', hue='model', style='model',
                 markers=True, lw=2)
    plt.legend(fontsize=12)
    plt.show()
