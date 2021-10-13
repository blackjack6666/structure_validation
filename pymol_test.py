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


def fasta_reader2(fasta_path:str):

    gene_protein_seq_dict = {}
    with open(fasta_path, 'r') as f_o:
        file_split = f_o.read().split('\n>')

    for each in file_split:
        first_line, seq = each.split('\n')[0], ''.join(each.split('\n')[1:])
        uniprot_id = first_line.split('|')[1]
        gene = first_line.split('GN=')[1].split(' ')[0] if 'GN=' in first_line else 'N/A'
        gene_protein_seq_dict[gene+'_'+uniprot_id] = seq
    return gene_protein_seq_dict


def fasta_reader3(fasta_path:str):

    protein_dict = {}
    with open(fasta_path, 'r') as f_o:
        file_split = f_o.read().split('\n>')

    for each in file_split:
        first_line, seq = each.split('\n')[0], ''.join(each.split('\n')[1:])
        uniprot_id = first_line.split('|')[1]
        gene = first_line.split('GN=')[1].split(' ')[0] if 'GN=' in first_line else 'N/A'
        des = ' '.join(first_line.split(' ')[1:]).split(' OS=')[0]
        protein_dict[uniprot_id] = (seq,gene,des)
    return protein_dict

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

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation

    for pep, new_pep, PTM in zip(peptide_list, new_pep_list,PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:

            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) -1
            # print (start_pos,end_pos,new_pep)
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
        matched_pep = tp[2]  # without ptm site
        zero_line[tp[0]:tp[1]+1]+=len(peptide_psm_dict[matched_pep])
        # print (tp[0],tp[1],matched_pep, 'aho')
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


def freq_ptm_index_gen_batch_v2(psm_list, protein_dict, regex_dict=None):
    """

    :param psm_list:
    :param protein_dict:
    :param regex_dict: {regex:HEX color}
    :return:
    """

    from collections import Counter
    from collections import defaultdict
    from heapq import heappush, heappop
    import time

    id_freq_array_dict = {}
    ptm_site_counting = defaultdict(int)
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    h = []
    regex_pat = '\w{1}\[\d+\.?\d+\]' # universal ptm pattern
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    for each in psm_list:

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
        id_freq_array_dict[id_list[i]] = zero_line_slice


        # if ptm_index_line_dict:
        #     id_ptm_idx_dict[id_list[i]]= {ptm:np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i]+1:separtor_pos_array[i+1]])[0]+1
        #                                   for ptm in ptm_index_line_dict}

        if ptm_index_line_dict:
            id_ptm_idx_dict[id_list[i]] = {ptm: np.array(
                np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[
                    0] + 1).tolist()
                                           for ptm in ptm_index_line_dict}
    print (time.time()-time_start)


    return id_freq_array_dict, id_ptm_idx_dict, [heappop(h) for i in range(len(h))][::-1]


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
    # pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
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
    # stick cysteines


    pymol.cmd.select('cysteines', 'resn cys')
    pymol.cmd.show('sticks', 'cysteines')
    pymol.cmd.set('stick_color', 'blue')
    pymol.cmd.set('stick_radius', 0.7)


    if png_sava_path:
        pymol.cmd.png(png_sava_path)

    print (f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    # dump_rep(pdb_name,base_path)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # pymol.cmd.save('new_file.pse')
    # Get out!
    # pymol.cmd.quit()


def show_cov_3d_v2(protein_id,
                   pdb_file,
                   id_freq_array_dict,
                   id_ptm_idx_dict,
                   regex_color_dict= None,
                   png_sava_path=None,
                   base_path=None):
    """

    :param protein_id:
    :param pdb_file:
    :param id_freq_array_dict: returned by freq_ptm_index_gen_batch_v2
    :param id_ptm_idx_dict: returned by freq_ptm_index_gen_batch_v2
    :param png_sava_path:
    :param base_path: html output base path
    :param regex_color_dict {regex: RGB_list}
    :return:
    """
    time_start = time.time()
    frequency_array = id_freq_array_dict[protein_id]
    if id_ptm_idx_dict != {}:
        ptm_nonzero_idx_dict = id_ptm_idx_dict[protein_id]

    else:
        ptm_nonzero_idx_dict = None


    pdb_name = os.path.split(pdb_file)[1]
    print(pdb_name)
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


    # # set customized color
    # if regex_color_dict:
    #     for i in regex_color_dict:
    #         pymol.cmd.set_color(i,regex_color_dict[i])
    #
    #
    # print (ptm_nonzero_idx_dict)
    # max_freq = np.max(frequency_array)
    # for i in range(len(frequency_array)): # iterate over each residue position
    #     # check if there is ptm
    #     ptm = False
    #     if ptm_nonzero_idx_dict:
    #
    #         for ptm_regex in ptm_nonzero_idx_dict:
    #             # if index has ptm
    #             if i in ptm_nonzero_idx_dict[ptm_regex]:
    #                 ptm = True
    #                 pymol.cmd.color(ptm_regex, 'resi %i' % (i + 1))
    #
    #     if ptm == False:  # if no ptm found
    #         if frequency_array[i] == 0:
    #             pymol.cmd.color('grey', 'resi %i' % (i + 1))
    #         elif 1 <= frequency_array[i] < 0.2 * max_freq:
    #             pymol.cmd.color('paleyellow', 'resi %i' % (i + 1))
    #         elif 0.2 * max_freq <= frequency_array[i] < 0.4 * max_freq:
    #             pymol.cmd.color('tv_yellow', 'resi %i' % (i + 1))
    #         elif 0.4 * max_freq <= frequency_array[i] < 0.6 * max_freq:
    #             pymol.cmd.color('yelloworange', 'resi %i' % (i + 1))
    #         elif 0.6 * max_freq <= frequency_array[i] < 0.8 * max_freq:
    #             pymol.cmd.color('tv_orange', 'resi %i' % (i + 1))
    #         else:
    #             pymol.cmd.color('sand', 'resi %i' % (i + 1))
    #     else: # ptm color assigned, move on to the next residue
    #         continue

    if png_sava_path:
        pymol.cmd.png(png_sava_path)

    print (f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    dump_rep_color_from_array(pdb_name,frequency_array,ptm_nonzero_idx_dict,regex_color_dict,base_path)
    # dump_rep(pdb_name,base_path)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # Get out!
    pymol.cmd.quit()


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
            # stick cysteines
            # pymol.cmd.select('cysteines','resn cys')
            # pymol.cmd.show('sticks','cysteines')

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
        pymol.cmd.quit()


def pdb_file_reader(pdb_file_list:list):
    """
    reads a pdb file into protein sequence
    :param pdb_file:
    :return:
    """
    import re
    from params import aa_dict,aa_reg_str
    import os

    pdb_protein_seq_dict = {}

    for pdb_file in pdb_file_list:
        with open(pdb_file,'r') as f_o:
            f_split = f_o.read().split('\nATOM')[1:]
            pos_aa_list = [(int(re.search('\d+(?=\s+[+-]?\d+\.)',each).group()),
                           re.search(aa_reg_str,each).group(0)) for each in f_split]
            protein_seq = ''
            for i in range(len(pos_aa_list)-1):
                if pos_aa_list[i+1][0] == pos_aa_list[i][0]:
                    continue
                else:
                    protein_seq += aa_dict[pos_aa_list[i][1]]

        # add last aa
        protein_seq+=aa_dict[pos_aa_list[-1][1]]
        pdb_protein_seq_dict[os.path.split(pdb_file)[-1]] = protein_seq
        print (protein_seq)
    return pdb_protein_seq_dict


if __name__=='__main__':

    import time
    import pandas as pd

    pdb_file_base_path = 'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/'

    peptide_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_240D/peptide.tsv'
    time_point_rep = ['120D','240D','1080D']
    psm_tsv_list = ['D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'/psm.tsv' for each in time_point_rep]
    print (f'{len(psm_tsv_list)} psm files to read...')
    fasta_file = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    peptide_list = peptide_counting(peptide_tsv)
    protein_dict = fasta_reader(fasta_file)
    psm_list = [psm for file in psm_tsv_list for psm in modified_peptide_from_psm(file)]

    # protein_list = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F_average.xlsx',index_col=0).index.tolist()
    protein_list = ['Q01149','P10853']
    print (protein_dict['Q01149'][385])


    protein_dict_sub = {prot:protein_dict[prot] for prot in protein_list}
    base_path = 'C:/tools/pymol-exporter-0.01/pymol_exporter/'  # JS folder path includes JS files required in html

    pdb_path_list = [file for each in protein_list for file in glob(pdb_file_base_path+'*'+each+'*.pdb')]
    # show_3d_batch(psm_list,protein_dict_sub,pdb_file_base_path,glmol_basepath=None)

    # show_cov_3d(psm_list,protein_dict_sub['P10853'],'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/AF-P10853-F1-model_v1.pdb',base_path='D:/data/alphafold_pdb/')


    # psm_list = ['FQSSAVM[147]ALQEACEAYLVGLFEDTNLCAIHAK']
    id_freq_array_dict, id_ptm_idx_dict, heap_list = freq_ptm_index_gen_batch_v2(psm_list,protein_dict,
                                                                                 regex_dict={'P\[113\]':[0,255,255], 'K\[144\]':[255,1,1]})
    # print (id_ptm_idx_dict['P68433'])
    # print (heap_list[0])
    show_cov_3d_v2('Q8TER0','D:/data/alphafold_pdb/UP000000589_10090_MOUSE/AF-Q8TER0-F1-model_v1.pdb',
                   id_freq_array_dict,id_ptm_idx_dict,regex_color_dict={'P\[113\]':[0,255,255], 'K\[144\]':[255,1,1]})
    

    # pdb_path = 'C:/Users/gao lab computer/Downloads/MS_SNED1_Suppl_File_S1.pdb'
    # pdb_path2 = 'D:/data/alphafold_pdb/AF-P11276-F1-model_v1.pdb'
    # pdb_file_reader([pdb_path])
    # fasta_file = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
    # protein_dict = fasta_reader3(fasta_file)
    # print (protein_dict['Q91W20'])

