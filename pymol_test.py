import re
import pymol
import numpy as np
from collections import defaultdict
import sys
from pymol2glmol import *
from glob import glob
import ahocorasick
import commons


def automaton_trie(peptide_list):
    A = ahocorasick.Automaton()
    for idx, peptide in enumerate(peptide_list):
        A.add_word(peptide, (idx, peptide))
    A.make_automaton()
    return A

def automaton_matching(A, seq_line):
    result = []
    for end_idx, (insert_order, original_value) in A.iter(seq_line):
        start_idx = end_idx - len(original_value) + 1
        result.append((start_idx, end_idx, original_value))
        assert seq_line[start_idx:start_idx+len(original_value)] == original_value
    return result


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def peptide_counting(peptide_tsv_file):

    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string,regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map single protein seq
    :param peptide_list:
    :param protein_seq_string:
    :param regex_pat:
    :return:
    """

    freq_array = np.zeros(len(protein_seq_string))
    PTM_sites_counting = defaultdict(int)
    PTM_loc_list = []
    #print peptide_list

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation
    for pep, new_pep, PTM in zip(peptide_list, new_pep_list,PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:
            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) -1
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:

                    PTM_index = pep.find(ele)
                   #PTM_site = pep[PTM_index] # single amino acid
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos+PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)
    return freq_array, PTM_loc_list, PTM_sites_counting


def freq_ptm_index_gen_batch(psm_list, protein_dict, regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map large psm list on whole proteome
    :param psm_list: psm list with ptms
    :param protein_dict: protein sequence dictionary
    :param regex_pat: regex pattern
    :return:
    """
    from collections import Counter
    from collections import defaultdict
    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}

    psm_count_dict = Counter(psm_list)  # count psm number for same sequence
    peptide_psm_dict = defaultdict(list) # append all psm into a dictionary
    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    ptm_index_line = commons.zero_line_for_seq(seq_line)
    separtor_pos_array = commons.separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]),seq_line)

    for tp in aho_result:
        matched_pep = tp[2]
        zero_line[tp[0]:tp[1]+1]+=psm_count_dict[matched_pep]
        for psm in peptide_psm_dict[matched_pep]:
            psm_mod = re.findall(regex_pat,psm)
            if psm_mod: # if psm has mod
                for ele in psm_mod:
                    ptm_idx = psm.find(ele)
                    ptm_index_line[tp[0]+ptm_idx]+=1

    for i in range(len(separtor_pos_array)-1):

        id_freq_array_dict[id_list[i]] = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]].tolist()
        id_ptm_idx_dict[id_list[i]] = np.nonzero(ptm_index_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1

    return id_freq_array_dict,id_ptm_idx_dict

def modified_peptide_from_psm(psm_path):
    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]',line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list



def show_cov_3d(peptide_list, protein_seq, pdb_file, png_sava_path=None, base_path=None):
    """
    show 3d coverage map based on peptide list and a protein_seq
    :param peptide_list:
    :param protein_seq:
    :param pdb_file:
    :param png_sava_path:
    :return:
    """
    import time
    time_start = time.time()
    freq_array, ptm_loc_list, PTM_sites_counting = freq_array_and_PTM_index_generator(peptide_list,protein_seq)

    print (ptm_loc_list)
    print (f'ptm sites counting: {PTM_sites_counting}')
    # open pdb file with pymol
    pdb_name = os.path.split(pdb_file)[1]
    print (pdb_name)
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable()
    print(pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('black')

    # highlight covered region
    max_freq = np.max(freq_array)
    for i, j in enumerate(freq_array):
        if i not in ptm_loc_list:
            if freq_array[i] == 0:
                pymol.cmd.color('grey', 'resi %i' % (i + 1))
            elif 1 <= freq_array[i] < 0.2 * max_freq:
                pymol.cmd.color('paleyellow', 'resi %i' % (i + 1))
            elif 0.2 * max_freq <= freq_array[i] < 0.4 * max_freq:
                pymol.cmd.color('tv_yellow', 'resi %i' % (i + 1))
            elif 0.4 * max_freq <= freq_array[i] < 0.6 * max_freq:
                pymol.cmd.color('yelloworange', 'resi %i' % (i + 1))
            elif 0.6 * max_freq <= freq_array[i] < 0.8 * max_freq:
                pymol.cmd.color('tv_orange', 'resi %i' % (i + 1))
            else:
                pymol.cmd.color('sand', 'resi %i' % (i + 1))
        else:
            pymol.cmd.color('red', 'resi %i' % (i + 1))

    if png_sava_path:
        pymol.cmd.png(png_sava_path)

    print (f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    dump_rep(pdb_name,base_path)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # pymol.cmd.save('new_file.pse')
    # Get out!
    # pymol.cmd.quit()


def show_3d_batch(psm_list, protein_dict, pdb_base_path, glmol_basepath=None):
    """
    output 3d html glmol in batch
    :param psm_list: psm list with modification M[#]
    :param protein_dict: protein dictionary
    :param pdb_file:
    :param base_path:
    :return:
    """
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()

    id_freq_array_dict, id_ptm_idx_dict = freq_ptm_index_gen_batch(psm_list,protein_dict)

    for id in id_freq_array_dict:
        time_start = time.time()
        pdb_path = pdb_base_path+'AF-'+id+'-F1-model_v1.pdb'
        if os.path.exists(pdb_path): # if local pdb db has this id entry
            pdb_name = os.path.split(pdb_path)[1]
            print(pdb_name)
            pymol.cmd.load(pdb_path, pdb_name)
            pymol.cmd.disable("all")
            pymol.cmd.enable()
            print(pymol.cmd.get_names())
            pymol.cmd.hide('all')
            pymol.cmd.show('cartoon')
            pymol.cmd.set('ray_opaque_background', 0)
            pymol.cmd.bg_color('black')

            # highlight covered region
            freq_array = id_freq_array_dict[id]
            ptm_loc_list = id_ptm_idx_dict[id]
            max_freq = np.max(freq_array)
            for i, j in enumerate(freq_array):
                if i not in ptm_loc_list:
                    if freq_array[i] == 0:
                        pymol.cmd.color('grey', 'resi %i' % (i + 1))
                    elif 1 <= freq_array[i] < 0.2 * max_freq:
                        pymol.cmd.color('paleyellow', 'resi %i' % (i + 1))
                    elif 0.2 * max_freq <= freq_array[i] < 0.4 * max_freq:
                        pymol.cmd.color('tv_yellow', 'resi %i' % (i + 1))
                    elif 0.4 * max_freq <= freq_array[i] < 0.6 * max_freq:
                        pymol.cmd.color('yelloworange', 'resi %i' % (i + 1))
                    elif 0.6 * max_freq <= freq_array[i] < 0.8 * max_freq:
                        pymol.cmd.color('tv_orange', 'resi %i' % (i + 1))
                    else:
                        pymol.cmd.color('sand', 'resi %i' % (i + 1))
                else:
                    pymol.cmd.color('red', 'resi %i' % (i + 1))


            # pymol2glmol, convert pdb to pse and visualize through html
            dump_rep(pdb_name, glmol_basepath)
            print(f'time used for mapping: {pdb_name, time.time() - time_start}')
        else:
            print (f'{id} not exist in pdb database')
            continue


if __name__=='__main__':
    import time
    import pandas as pd

    pdb_file_base_path = 'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/'

    peptide_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_240D/peptide.tsv'
    psm_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_240D/psm.tsv'
    fasta_file = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    peptide_list = peptide_counting(peptide_tsv)
    protein_dict = fasta_reader(fasta_file)
    psm_list = modified_peptide_from_psm(psm_tsv)

    protein_list = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F_average.xlsx',index_col=0).index.tolist()
    protein_dict_sub = {prot:protein_dict[prot] for prot in protein_list}
    base_path = 'C:/tools/pymol-exporter-0.01/pymol_exporter/'  # JS folder path required for html

    pdb_path_list = [file for each in protein_list for file in glob(pdb_file_base_path+'*'+each+'*.pdb')]
    show_3d_batch(psm_list,protein_dict_sub,pdb_file_base_path,base_path)
    # for protein_id, pdb in zip(protein_list,pdb_path_list):
    #
    #     show_cov_3d(psm_list,protein_dict[protein_id],pdb, base_path=base_path)

    # pymol.cmd.load('your_session.pse')
    # dump_rep('AF-Q5ZQU0-F1-model_v1',base_path='C:/tools/pymol-exporter-0.01/pymol_exporter')