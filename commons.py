from collections import defaultdict
import numpy as np
import glob
import pickle as ppp
import multiprocessing as mp
import time

import re


def extract_UNID_and_seq(protein_dict):
    UNID_list = [key for key in protein_dict.keys()]
    seq_list = [value for value in protein_dict.values()]
    return UNID_list, seq_list

def read_seq_length_into_nparray(seq_list):
    seq_length_array = np.array([])
    for each_seq in seq_list:
        seq_length_array = np.append(seq_length_array, len(each_seq))
    return seq_length_array

def creat_total_seq_line(seq_list, sep=None):
    seq_line = '|'.join(seq_list) if sep == '|' else ''.join(seq_list)
    return seq_line

def zero_line_for_seq(seq_line):
    zero_line = np.zeros(len(seq_line),dtype=np.int32)
    return zero_line

# the following function create a dictionary that read each position in long sequence line as key, corresponding UNIPORT ID as value.
def read_position_ID_into_dict(UNID_list, seq_list, zero_line):
    m = 0
    j = 0
    seq_line_ID_dict = dict()
    for i in range(len(zero_line)):
        if j < len(seq_list[m]):
            seq_line_ID_dict[i] = UNID_list[m]
            j += 1
        else:
            j = 0

            m += 1
    return seq_line_ID_dict


def creat_ID_pep_dict(aho_result, pos_ID_dict):
    ID_pep_dict = defaultdict(set)
    for i in aho_result:
        ID_pep_dict[pos_ID_dict[i[0]]].add(i[2])
    return ID_pep_dict


def creat_pep_ID_dict(aho_result, pos_ID_dict):
    pep_ID_dict = defaultdict(set)
    for i in aho_result:
        pep_ID_dict[i[2]].add(pos_ID_dict[i[0]])
    return pep_ID_dict


def create_unique_id_peptide_dict(pep_id_dict):
    """
    get a dictionary with unique peptides for each protein
    :param pep_id_dict:
    :return:
    """
    unique_id_peptide_dict = defaultdict(set)
    unique_id_peptide_count_dict = defaultdict(int)
    unique_pep_id_dict = {pep:prot for pep in pep_id_dict for prot in pep_id_dict[pep]
                          if len(pep_id_dict[pep])==1}
    for pep in unique_pep_id_dict:
        unique_id_peptide_dict[unique_pep_id_dict[pep]].add(pep)

    for id in unique_id_peptide_dict:
        unique_id_peptide_count_dict[id]=len(unique_id_peptide_dict[id])

    return unique_id_peptide_dict, unique_id_peptide_count_dict



def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def find_all(a_str, sub):  # a generator
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start +=  len(sub)  # use start += 1 to find overlapping matches

def find_pep_start(param):
    seq_line, peptide_list=param
    res_dict={}
    for peptide in peptide_list:
        res_dict[peptide]=[m for m in find_all(seq_line,peptide)]
    return res_dict

# the following function returns a dictionary with each peptide as key and corresponding start position list as value.
def start_end_pos_dict(res_dicts):
    start_pos_dict = {}
    end_pos_dict = {}
    for res_dict in res_dicts:
        for peptide in res_dict:
            start_pos_dict[peptide] = res_dict[peptide]
            end_pos_dict[peptide] = [i + len(peptide) - 1 for i in res_dict[peptide]]
    return start_pos_dict, end_pos_dict

# the following function returns a number_line after matching peptides to long_seq_line.
def adding_numbers_to_zero_line(zero_line, start_dict, end_dict):
    for peptide in start_dict:
        start_list_for_each = start_dict[peptide]
        end_list_for_each = end_dict[peptide]
        for i, j in zip(start_list_for_each, end_list_for_each):
            zero_line[i:j] += 1
    return zero_line

def separator_pos(seq_line):
    sep_pos_array = np.array([m.start() for m in re.finditer('\|', seq_line)],dtype=np.int32)
    sep_pos_array = np.insert(sep_pos_array, 0, 0)
    sep_pos_array = np.append(sep_pos_array, len(seq_line))
    return sep_pos_array


def get_unique_peptide(list_of_peptsv:list):
    """
    from pep tsv file only get unique peptides compared with previous ones, e.g. in 4 hour sample, filter out peptides
    in 1h,2h and only retain peptides uniquely identified in 4h
    :param list_of_peptide:
    :return:
    """
    from pymol_test import peptide_counting
    unique_peptide_dict = {}
    peptide_list = []
    for idx, val in enumerate(list_of_peptsv):
        if '\\' in val:
            file_name = val.split('\\')[-2]
        else:
            file_name = val.split("/")[-2]
        print (file_name)
        unique_peptide_list = [each for each in peptide_counting(val) if each not in peptide_list]

        peptide_list += unique_peptide_list

        unique_peptide_dict[file_name] = unique_peptide_list

    return unique_peptide_dict


def miss_cleavage_identify(peptide_list,regex_pattern:dict):
    """
    calculate the ratio of peptides that has least one missed cleavage to all peptides
    :param peptide_list:
    :param regex_pattern: regex pattern for mapping enzyme specificity. Default is trypsin
    :return:
    """
    from re import findall
    k_miss_sum,r_miss_sum = 0,0
    for each in peptide_list:
        k_miss_sum += len(findall(regex_pattern['K'],each))
        r_miss_sum += len(findall(regex_pattern['R'],each))
        # print (each,findall(regex_pattern['K'],each))
    return k_miss_sum,r_miss_sum


def dta_charge_reader(file_location_path):
    """
    -----
    read peptide seq, PSM, protein ID from dta files
    -----
    :param file_location_path: could be str folder path or dta file path or list of dta files or folder path
    :return: peptide list, PSM dictionary showing number of PSM, protein ID list
    """
    from glob import glob
    from collections import defaultdict

    dta_files = []

    # single dta read
    if isinstance(file_location_path,str) and \
            (file_location_path.endswith('.dta') or file_location_path.endswith('.txt')):
        dta_files = [file_location_path]
        print ('reading single dta file')

    # dta folder read
    elif isinstance(file_location_path,str) and not file_location_path.endswith('.dta'):
        dta_files = glob(file_location_path + '*.dta')
        print ('reading single folder')

    # multiple dta files read or multiple dta folder read
    elif isinstance(file_location_path,list):
        for each_dta in file_location_path:
            # when parameter is dta file

            if each_dta.endswith('.dta'):
               dta_files.append(each_dta)

            # when parameter is a folder path
            else:
                dta_files += glob(each_dta+'*.dta')
    else:
        raise ValueError('parameter should be string folder path or dta file path or list of dta files or folder paths')

    # exclude wash and hela files
    clean_dta_files = []
    for each in dta_files:
        wash_hela = 0
        for word in ['wash', 'Wash', 'WASH', 'Hela', 'hela', 'HELA']:
            if word in each:
                wash_hela += 1
                break
        if wash_hela == 0:
            clean_dta_files.append(each)

    print (clean_dta_files)

    # read info.

    seq_charge_dict = defaultdict(list)
    for dta_file in clean_dta_files:
        with open(dta_file, 'r') as file_open:
            for i in range(38):
                next(file_open)
            Reverse_start = 0
            for line in file_open:
                line_split = line.split('\t')
                # print (len(line_split))
                if line.startswith('Reverse_') or line.startswith('Rev_'):
                    Reverse_start = 1
                elif line.startswith('sp') or line.startswith('tr'):
                    Reverse_start = 0

                elif len(line_split) == 14 and Reverse_start == 0:
                    # print (line_split[-1].split('.'))
                    pep_seq = line_split[-1].split('.')[1]
                    # print (pep_seq)
                    charge = int(line_split[1].split('.')[-1])
                    seq_charge_dict[pep_seq].append(charge)

    seq_charge_dict = {each:list(set(seq_charge_dict[each])) for each in seq_charge_dict}
    return seq_charge_dict

if __name__ == "__main__":
    """
    filename = 'C:/uic/lab/data/xinhao_data1/uniprot-proteome_UP000005640.fasta'
    test_filename = 'C:/uic/lab/data/TEST/test_fasta.txt'
    protein_dict = read_fasta_into_dict(filename)[0]
    #print protein_dict, protein_dict.values()
    uniprot_ID_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list)
    zero_line = zero_line_for_seq(seq_line)
    seq_line_ID_dict = read_position_ID_into_dict(uniprot_ID_list, seq_list, zero_line)
    #ppp.dump(seq_line_ID_dict, open('ID_position_dict.P'), protocol=-1) # read it into a pickle file
    path = 'C:/uic/lab/data/xinhao_data1/'
    test_path = 'C:/uic/lab/data/TEST/'
    dtafiles = glob.glob(path+'*.dta')
    start = time.clock()
    peptide_list = read_peptide(dtafiles)
    peptide_list_set = set(peptide_list)
    peptide_list_unique = list(peptide_list_set)
    print (time.clock()-start, len(peptide_list), len(peptide_list_unique))

    #chunk_list = chunks(peptide_list, 10)
    chunk_list = chunks(peptide_list_unique, 10)
    parm_list = [(seq_line, p) for p in chunk_list]
    start = time.clock()
    pool = mp.Pool(processes = mp.cpu_count()-4)
    res_dicts = pool.map(find_pep_start, parm_list)
    pool.close()
    pool.join()

    start_pos_dict, end_pos_dict = start_end_pos_dict(res_dicts)
    zero_line = adding_numbers_to_zero_line(zero_line, start_pos_dict, end_pos_dict)
    print (time.clock()-start)

    sep_pos_array = separator_pos(seq_line)
    # trie implementation
    '''
    zero_line_trie = in_trie(make_trie(peptide_list_unique), seq_line)
    total_seq_len_trie = len(zero_line_trie)-len(sep_pos_array) + 2
    total_non_zero_trie = np.count_nonzero(zero_line_trie)
    overall_percentage_trie = float(total_non_zero_trie)/total_seq_len_trie*100
    print total_non_zero_trie
    zero_line_naive = ppp.load(open('zero_line.p', 'rb'))
    print np.count_nonzero(zero_line_naive)
    '''
    print (time.clock()-start)
    total_seq_len = len(zero_line) - len(sep_pos_array) + 2  # artifically added two separator positions into sep_pos_array so need to plus 2
    total_non_zero = np.count_nonzero(zero_line)
    overall_percentage = float(total_non_zero) / total_seq_len * 100
    print (overall_percentage)
    #print len(uniprot_ID_list), len(sep_pos_array)
    ppp.dump(zero_line, open('zero_line.p', 'wb'), protocol=-1)
    ppp.dump(sep_pos_array, open('separator_pos.p', 'wb'), protocol=-1)
    ppp.dump(uniprot_ID_list, open('uniprotID.p', 'wb'), protocol=-1)
        #ppp.dump(uniprot_ID_list, pf)

    #print protein_dict, seq_line, peptide_list, sep_pos_array, zero_line_trie
    """
    # import os
    # base_path = 'D:/data/native_protein_digestion/11052021/search_result/'
    # peptide_tsv_list = [base_path + folder+'/peptide.tsv' for folder in os.listdir(base_path)
    #            if os.path.isdir(os.path.join(base_path, folder))]
    # print (peptide_tsv_list)
    # peptide_dict = get_unique_peptide(peptide_tsv_list)
    # print ((len(peptide_dict['0005min'])))

    from glob import glob
    import imageio
    filenames = glob('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/*unique.png')
    print (filenames)
    images = []
    with imageio.get_writer('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/Q96I99.gif',
                            mode='I',duration=1,fps=30) as writer:
        for file in filenames:
            image = imageio.imread(file)
            writer.append_data(image)

