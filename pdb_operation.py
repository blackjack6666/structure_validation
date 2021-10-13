### some functions to operate on PDB file

import re
from collections import defaultdict
from pymol_test import *
import numpy as np
import time


def pdb_file_reader(pdb_file):
    """
    read pdb file and map xyz coordinates of each residue
    :param pdb_file:
    :return:
    """
    with open(pdb_file,'r') as f_o:
        file_split = f_o.read().split('\nATOM')[1:]

    residue_atom_xyz = defaultdict(list)

    # append xyz coords of each atom to residue position
    for line in file_split:
        residue_atom_xyz[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append([float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)])
    # print (residue_atom_xyz)
    return residue_atom_xyz


def residue_distance(residue_atom_xyz):
    """
    calculate the distance for each residue to center of 3d structure
    :param residue_atom_xyz:
    :return:
    """
    zero_point = np.array([0,0,0])
    residue_distance_dict = {}
    for each_pos in residue_atom_xyz:
        distance = 0
        for each_atom in residue_atom_xyz[each_pos]:
            # print (each_atom)
            distance += np.linalg.norm(np.array(each_atom)-zero_point)
        distance = distance/len(residue_atom_xyz[each_pos])
        residue_distance_dict[each_pos] = distance
    return residue_distance_dict


def cov_distance(freq_array,residue_dist_dict):
    """
    calculate averge distance of covered region in one protein
    :param freq_array:
    :param residue_dist_dict:
    :return:
    """
    ave_dist = 0
    non_zero_index = np.nonzero(freq_array)[0]
    num_nonzeros = np.count_nonzero(freq_array)
    for i in non_zero_index:
        ave_dist += residue_dist_dict[i+1]

    if num_nonzeros == 0:
        return 100
    else:
        return ave_dist/num_nonzeros

if __name__ == '__main__':
    import time

    time_point_rep = ['1h','2h','4h','18h']
    psm_tsv_list = ['D:/data/native_protein_digestion/' + each + '_1_native/psm.tsv' for each in time_point_rep]
    print(f'{len(psm_tsv_list)} psm files to read...')
    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
    # peptide_list = peptide_counting(peptide_tsv)
    protein_dict = fasta_reader(fasta_file)
    psm_list = [psm for file in psm_tsv_list for psm in modified_peptide_from_psm(file)]

    pdb_file = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-P61604-F1-model_v1.pdb'
    time_start = time.time()
    residue_distance = residue_distance(pdb_file_reader(pdb_file))
    print (time.time()-time_start)
    print (residue_distance)
    freq_array = freq_ptm_index_gen_batch_v2(psm_list,protein_dict)[0]['P61604']

    print (cov_distance(freq_array,residue_distance))



