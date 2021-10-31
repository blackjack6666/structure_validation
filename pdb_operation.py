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


def find_centroid(residue_atom_xyz):
    """
    find geometric center of protein given xyz coordinates
    :param residue_atom_xyz:
    :return: use median as centroid
    """
    atom_coordinates = np.array([coord for each in residue_atom_xyz for coord in residue_atom_xyz[each]])

    return np.median(atom_coordinates,axis=0)


def residue_distance(residue_atom_xyz):
    """
    calculate the distance for each residue to center of 3d structure
    :param residue_atom_xyz:
    :return:
    """
    # zero_point = np.array([0,0,0])
    zero_point = find_centroid(residue_atom_xyz)
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

    num_nonzeros = np.count_nonzero(freq_array)
    if num_nonzeros == 0:
        return None
    else:
        ave_dist = 0
        non_zero_index = np.nonzero(freq_array)[0]
        for i in non_zero_index:
            ave_dist += residue_dist_dict[i+1]

        return ave_dist/num_nonzeros


def pdb2_3darray(pdb_file):
    """
    convert pdb coordinates into 3d numpy array
    :param pdb_file:
    :return:
    """
    protein_cood_dict = pdb_file_reader(pdb_file)
    normalized_protein_coord_dict = {}

    # average the coordinates of atoms
    for _ in protein_cood_dict:
        locs = np.mean(np.array(protein_cood_dict[_]),axis=0)
        protein_cood_dict[_] = locs

    new_coord_array = np.array([each for each in protein_cood_dict.values()])

    # get x, y, z range
    x_diff,y_diff,z_diff = [(i-j,j) for i, j in zip(np.max(new_coord_array, axis=0), np.min(new_coord_array,axis=0))]

    # normalize x y z coordinates
    for _ in protein_cood_dict:
        x,y,z = protein_cood_dict[_]
        new_x,new_y,new_z = int((x-x_diff[1])*10), int((y-y_diff[1])*10),int((z-z_diff[1])*10)
        print(new_x,new_y,new_z)
        normalized_protein_coord_dict[_-1] = (new_x,new_y,new_z)

    return normalized_protein_coord_dict, [int(x_diff[0]*10), int(y_diff[0]*10), int(z_diff[0]*10)]


def map_aa2_3darray(freq_array, normalized_protein_coord_dict, coord_range_list):
    """
    map aa
    :param freq_array:
    :param normalized_protein_coord_dict:
    :param coord_range_list: [x_range, y_range, z_range], returned by pdb2_3darray second return
    :return:
    """
    x,y,z = coord_range_list
    protein_3d = np.zeros((x+1,y+1,z+1),dtype=np.int8)
    non_zero_index = np.nonzero(freq_array)[0]

    for each in non_zero_index:
        # sequence index to 3d coordinates
        x,y,z = normalized_protein_coord_dict[each]

        protein_3d[x][y][z] +=1
    return protein_3d

if __name__ == '__main__':
    import time

    pdb_file = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-P61604-F1-model_v1.pdb'
    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
    protein_dict = fasta_reader(fasta_file)

    """
    time_point_rep = ['1h','2h','4h','18h']
    psm_tsv_list = ['D:/data/native_protein_digestion/' + each + '_1_native/psm.tsv' for each in time_point_rep]
    print(f'{len(psm_tsv_list)} psm files to read...')
    
    # peptide_list = peptide_counting(peptide_tsv)
    
    psm_list = [psm for file in psm_tsv_list for psm in modified_peptide_from_psm(file)]

    # time_start = time.time()
    # residue_distance = residue_distance(pdb_file_reader(pdb_file))
    # print (time.time()-time_start)
    # print (residue_distance)
    freq_array = freq_ptm_index_gen_batch_v2(psm_list,protein_dict)[0]['P61604']
    normalized_prot_dict, xyz_list = pdb2_3darray(pdb_file)
    print (xyz_list)
    mapped_3darray = map_aa2_3darray(freq_array,normalized_prot_dict,xyz_list)
    flat_array = mapped_3darray.flatten()
    print (np.count_nonzero(flat_array))
    #
    # print (cov_distance(freq_array,residue_distance))

    """

    ### calculate covered distance and write to excel
    def protein_tsv_reader(protein_tsv_file):
        with open(protein_tsv_file, 'r') as file_open:
            next(file_open)
            return [line.split("\t")[3] for line in file_open]

    protein_tsv = 'D:/data/native_protein_digestion/10282021/search_result_4miss/combined_protein.tsv'
    protein_list = protein_tsv_reader(protein_tsv)

    pdb_base = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
    time_point_rep = ['01h_h2o', '02h_h2o', '04h_h2o', '20h_h2o']
    psm_tsv_list = ['D:/data/native_protein_digestion/10282021/search_result_4miss/' + each + '/peptide.tsv' for each in time_point_rep]
    print(f'{len(psm_tsv_list)} psm files to read...')

    import pandas as pd
    df = pd.DataFrame(index=protein_list, columns=time_point_rep)

    for pep_tsv in psm_tsv_list:
        print (pep_tsv)
        peptide_list = peptide_counting(pep_tsv)
        freq_array_dict = freq_ptm_index_gen_batch_v2(peptide_list,protein_dict)[0]
        for prot in protein_list:
            pdb_file_path = pdb_base+'AF-'+prot+'-F1-model_v1.pdb'
            if os.path.exists(pdb_file_path):
                residue_dist_dict = residue_distance(pdb_file_reader(pdb_file_path))
                if len(residue_dist_dict) == len(protein_dict[prot]):

                    freq_array = freq_array_dict[prot]
                    cov_dist = cov_distance(freq_array,residue_dist_dict)
                    df.at[prot,pep_tsv.split('/')[-2]] = cov_dist
                else:
                    print ('%s protein len between pdb and fasta is not same' % prot)
            else:
                continue
    df.to_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist.xlsx')
