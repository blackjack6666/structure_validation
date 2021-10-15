import commons
from collections import defaultdict
import re
from pymol_test import automaton_matching, automaton_trie
import numpy as np
import time

def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_ptm_index_gen_batch_v2(psm_list_dict, protein_dict, regex_dict=None):
    """

    :param psm_list_dict: a dictionary with sample name as key and psm list/peptide list as value, ex. {1h:peplist1,
    2h:peplist2...}
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """

    from collections import Counter
    from collections import defaultdict
    from heapq import heappush, heappop

    sample_id_freq_dict = {}
    sample_id_ptm_idx_dict = {}
    sample_heap_dict = {}
    for sample in psm_list_dict:
        id_freq_array_dict = {}
        id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
        h = []
        regex_pat = '\w{1}\[\d+\.?\d+\]'  # universal ptm pattern
        peptide_psm_dict = defaultdict(list)

        for each in psm_list_dict[sample]:

            each_reg_sub = re.sub(regex_pat, my_replace, each)
            peptide_psm_dict[each_reg_sub].append(each)

        # aho mapping
        id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
        seq_line = commons.creat_total_seq_line(seq_list, sep="|")
        zero_line = commons.zero_line_for_seq(seq_line)
        ptm_index_line_dict = {each:commons.zero_line_for_seq(seq_line) for each in regex_dict} if regex_dict else False
        separtor_pos_array = commons.separator_pos(seq_line)

        aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)
        for tp in aho_result:
            matched_pep = tp[2]  # without ptm site
            zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
            # ptm assign might need to be optimized
            if ptm_index_line_dict:  # if ptm enabled
                for psm in peptide_psm_dict[matched_pep]:
                    for ptm in regex_dict:
                        # only keep one ptm in psm if there are multiple for correct index finding
                        new_psm = re.sub('(?!'+ptm[1:]+')\[\d+\.?\d+\]','',psm)

                        ptm_mod = re.findall(ptm, new_psm)
                        print(ptm_mod)
                        if ptm_mod:
                            for ele in ptm_mod:
                                print (new_psm)
                                ptm_idx = new_psm.find(ele)
                                print(matched_pep, tp[0], ele, ptm_idx)
                                ptm_index_line_dict[ptm][tp[0] + ptm_idx] += 1

        time_start = time.time()
        for i in range(len(separtor_pos_array)-1):
            zero_line_slice = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]]
            percentage_cov = np.count_nonzero(zero_line_slice)/len(zero_line_slice)*100
            heappush(h,(percentage_cov,id_list[i],zero_line_slice))
            id_freq_array_dict[id_list[i]] = zero_line_slice.tolist()


            # if ptm_index_line_dict:
            #     id_ptm_idx_dict[id_list[i]]= {ptm:np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1
            #                                   for ptm in ptm_index_line_dict}

            if ptm_index_line_dict:
                id_ptm_idx_dict[id_list[i]] = {ptm: np.array(
                    np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[
                        0] + 1).tolist()
                                               for ptm in ptm_index_line_dict}
        sample_id_freq_dict[sample] = id_freq_array_dict
        sample_id_ptm_idx_dict[sample] = id_ptm_idx_dict
        sample_heap_dict[sample] = [heappop(h) for i in range(len(h))][::-1]

    return sample_id_freq_dict,sample_id_ptm_idx_dict,sample_heap_dict


def show_3d_multiple(protein_id,
                     pdb_file,
                     sample_id_freq_dict,
                     sample_id_ptm_idx_dict,
                     regex_dict=None):
    """
    return pymol instance that has multiple structures
    :param protein_id:
    :param pdb_file:
    :param sample_id_freq_dict: {1h:{id1:[],id2:[]...},2h:{id1:[],id2:[]...}...}
    :param sample_id_ptm_idx_dict:
    :param regex_dict:
    :return:
    """


def pymol_imple(pdb_file_list:list):

    import pymol
    import os

    for pdb_file in pdb_file_list:
        pdb_name = os.path.split(pdb_file)[1]
        print(pdb_name)
        pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
        pymol.finish_launching()
        pymol.cmd.load(pdb_file, pdb_name)
