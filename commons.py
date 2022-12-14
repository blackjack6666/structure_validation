from collections import defaultdict
import numpy as np
import glob
import pickle as ppp
import multiprocessing as mp
import time

import re
aa_mass_table = {'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943,
                 'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464,
                 'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963,
                 'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028,
                 'T': 101.047679, 'U': 150.95363, 'W': 186.079313, 'Y': 163.06332,
                 'V': 99.068414, 'B': 114.53494, 'Z': 128.55059, 'X':110,
                }

expasy_rules = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3': r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                    r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                    r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV])',
    'thrombin': r'((?<=G)R(?=G))|'
                r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
}

def protein_mass(protein_seq):
    mass = 0
    for aa in protein_seq:
        mass += aa_mass_table[aa]
    mass += 18.01528 # add the weight of water

    return mass

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


def zero_line_for_seq(seq_line, datatype=np.int32):
    zero_line = np.zeros(len(seq_line), dtype=datatype)
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


def get_aggre_peptide(list_of_peptsv: list):
    """
    aggregate peptides at each time point
    :param list_of_peptsv:
    :return:
    """
    from pymol_test import peptide_counting
    aggre_peptide_dict = {}
    peptide_list = []
    for idx, val in enumerate(list_of_peptsv):
        if '\\' in val:
            file_name = val.split('\\')[-2]
        else:
            file_name = val.split("/")[-2]
        print(file_name)
        unique_peptide_list = [each for each in peptide_counting(val) if each not in peptide_list]

        peptide_list += unique_peptide_list
        print(file_name, len(peptide_list))
        aggre_peptide_dict[file_name] = set(peptide_list)
    print([(each, len(aggre_peptide_dict[each])) for each in aggre_peptide_dict])
    return aggre_peptide_dict


def psm_reader(psm_path, fragpipe_ver=13.0):
    pep_spec_count_dict = defaultdict(int)
    ret_pep_dict = {}
    with open(psm_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pep_seq = line_split[2] if fragpipe_ver == 13.0 else line_split[1]
            # retention_time = float(line_split[5])/60  if fragpipe_ver==13.0 else float(line_split[4])/60 # in minute
            pep_spec_count_dict[pep_seq] += 1
            # ret_pep_dict[retention_time] = pep_seq
    return pep_spec_count_dict

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


def dta_reader(dta_file_list):
    peptide_dict = defaultdict(list)
    for dta_file in dta_file_list:
        with open(dta_file, 'r') as file_open:
            for i in range(31):
                next(file_open)
            Reverse_start = 0
            for line in file_open:
                line_split = line.split('\t')

                if line.startswith('Reverse_') or line.startswith('Rev_'):
                    Reverse_start = 1
                elif line.startswith('sp') or line.startswith('tr'):
                    Reverse_start = 0
                    protein_id = line_split[0].split('|')[1]
                elif len(line_split) == 15 and Reverse_start == 0:
                    file = line_split[1].split('_')[0]
                    if file == 'GAPDH' or 'HGAPDH' or 'MGAPDH':
                        # print (line_split[-1])
                        pep_seq = line_split[-1].split('.')[1]
                        pep_mass = line_split[5]
                        if 'K' in pep_seq:
                            peptide_dict[file].append((pep_seq,pep_mass))
                        # elif len(re.findall('K',pep_seq))==2:
                        #     print(file,pep_seq,pep_mass)
    return peptide_dict


def dta_reader2(dta_file_list):
    peptide_mass_dict = defaultdict(list)
    for dta_file in dta_file_list:
        protein_id_list = []
        with open(dta_file, 'r') as file_open:
            for i in range(43):
                next(file_open)
            Reverse_start = 0
            pep_seq_switch = 0
            for line in file_open:
                line_split = line.split('\t')

                if line.startswith('Reverse_') or line.startswith('Rev_'):
                    Reverse_start = 1
                elif (line.startswith('sp') or line.startswith('tr')) and pep_seq_switch==0:
                    Reverse_start = 0
                    protein_id = line_split[0].split('|')[1]
                    # print (protein_id)
                    protein_id_list.append(protein_id)
                elif (line.startswith('sp') or line.startswith('tr')) and pep_seq_switch==1:
                    protein_id_list = []
                    Reverse_start = 0
                    protein_id = line_split[0].split('|')[1]
                    protein_id_list.append(protein_id)
                    pep_seq_switch = 0
                elif len(line_split) == 15 and Reverse_start == 0:  #  sometimes len is 14
                    pep_seq_switch = 1
                    pep_seq = line_split[-1].split('.')[1]
                    pep_mass = line_split[5]
                    file_name = line_split[1].split('_')[0]
                    for protein_id in protein_id_list:
                        peptide_mass_dict[protein_id].append((pep_seq,pep_mass,file_name))
    return peptide_mass_dict


def dta_result_analyze(peptide_mass_dict):
    """
    distinguish between light and heavy labeled peptide
    :param peptide_mass_dict: returned from dta_reader_2
    :return:
    """
    from collections import Counter
    prot_iso_lab_dict = {}
    for prt in peptide_mass_dict:
        isolabel_dict = defaultdict(set)
        pep_mass_lst = peptide_mass_dict[prt]
        for tp in pep_mass_lst:
            peptide_seq, observed_mass,file_name = tp

            # therotical_mass = protein_mass(peptide_seq)
            # num_of_k = Counter(peptide_seq)['K']
            # mass_diff = ((observed_mass-therotical_mass)-36)/num_of_k

            if file_name[0]=='H':  # probably means heavy
                isolabel_dict['heavy'].add(peptide_seq)
            else:
                isolabel_dict['light'].add(peptide_seq)
        prot_iso_lab_dict[prt] = isolabel_dict
    return prot_iso_lab_dict


def protein_tsv_reader(protein_tsv_file, protein_column=1):
    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        return [line.split("\t")[protein_column] for line in file_open]


def protein_info_from_fasta(fasta_path):
    """
    get protein name, gene name, entry name, and description
    :param fasta_path:
    :return:
    """
    info_dict = {}
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line.split('|')[1]
                cls = line.split('|')[0].split('>')[1]
                # print (protein_id)
                description = ' '.join(line.split('OS=')[0].split(' ')[1:])

                gene_name = line.split('GN=')[1].split(' ')[0].rstrip('\n') if 'GN=' in line else 'N/A'
                info_dict[protein_id] = (gene_name, description, cls)
    return info_dict
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

    # from glob import glob
    # import imageio
    # filenames = glob('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/*unique.png')
    # print (filenames)
    # images = []
    # with imageio.get_writer('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/Q96I99.gif',
    #                         mode='I',duration=1,fps=30) as writer:
    #     for file in filenames:
    #         image = imageio.imread(file)
    #         writer.append_data(image)

    dta_file = 'D:/data/native_protein_digestion/dimethylation/DTASELECT/030515_control_HL_PBS_chy_corr_2015_04_14_14_31470_DTASelect-filter.txt'
    file_list = glob.glob('D:/data/native_protein_digestion/dimethylation/DTASELECT/1_*DTASelect-filter.txt')+\
                glob.glob('D:/data/native_protein_digestion/dimethylation/DTASELECT/2_*DTASelect-filter.txt')

    peptide_mass_dict = dta_reader2([dta_file])

    prot_iso_label_dict = dta_result_analyze(peptide_mass_dict)['P84077']
    print (prot_iso_label_dict)