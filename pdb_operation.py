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


def find_centroid(residue_atom_xyz, centroid_method='mean'):
    """
    find geometric center of protein given xyz coordinates
    :param residue_atom_xyz:
    :return: use median as centroid
    """
    atom_coordinates = np.array([coord for each in residue_atom_xyz for coord in residue_atom_xyz[each]])

    if centroid_method == 'mean':
        return np.mean(atom_coordinates,axis=0)
    elif centroid_method == 'median':
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

        total_dist = sum([np.linalg.norm(np.array(each_atom)-zero_point)
                          for each_atom in residue_atom_xyz[each_pos]])

        average_dist = total_dist/len(residue_atom_xyz[each_pos])
        residue_distance_dict[each_pos] = average_dist
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


def cov_dist_normalize(freq_array,residue_dist_dict):
    """
    calculate the normalized distance of covered region in one protein, normalized_dist = dist/max(dist)
    :param freq_array:
    :param residue_dist_dict:
    :return:
    """
    num_nonzeros = np.count_nonzero(freq_array)
    max_dist = max([v for v in residue_dist_dict.values()])
    if num_nonzeros == 0:
        return None
    else:
        ave_dist = 0
        non_zero_index = np.nonzero(freq_array)[0]
        for i in non_zero_index:
            ave_dist += residue_dist_dict[i + 1]/max_dist*100

        return ave_dist / num_nonzeros


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
    psm_tsv_list = ['D:/data/native_protein_digestion/10282021/search_result_4miss/h20/' + each + '/peptide.tsv' for each in time_point_rep]
    print(f'{len(psm_tsv_list)} psm files to read...')
    """
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
    df.to_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist_centroid_mean.xlsx')
    """

    import matplotlib.pyplot as plt
    from statistics import mean
    import random
    # prots_tocheck = random.sample(protein_list,100)
    prots_tocheck = [each.split('\\')[-1].split('.png')[0] for each in glob('D:/data/native_protein_digestion/10282021/protein_centroid/*')]
    print (prots_tocheck)

    ### plot 3d and centroid
    """
    for prot in prots_tocheck:
        pdb_file_path = pdb_base+'AF-'+prot+'-F1-model_v1.pdb'

        if os.path.exists(pdb_file_path):
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            print (prot)
            residue_atom_xyz = pdb_file_reader(pdb_file_path)
            xyz = [each for v in residue_atom_xyz.values() for each in v]
            x,y,z = zip(*xyz)
            ax.scatter(x,y,z,marker='o',s=0.5)
            centroid = find_centroid(residue_atom_xyz)
            # plot centroid
            ax.scatter([centroid[0]],[centroid[1]],[centroid[2]], marker='o', s=8,color='r')
            ax.text2D(0.05, 0.95, "centroid coordinates: %.2f,%.2f,%.2f" % (centroid[0] ,centroid[1],centroid[2]), transform=ax.transAxes)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            plt.savefig('D:/data/native_protein_digestion/10282021/protein_centroid_median/%s.png' % prot)
    """
    ### plot residue distance distribution
    bot,med,top,all = [],[],[],[]
    for prot in prots_tocheck:
        pdb_file_path = pdb_base + 'AF-' + prot + '-F1-model_v1.pdb'

        if os.path.exists(pdb_file_path):
            distance_array = sorted([v for v in residue_distance(pdb_file_reader(pdb_file_path)).values()])
            block = int(len(distance_array)*0.33)
            top_33,med_33,bot_33 = distance_array[-block:],distance_array[block:-block],distance_array[:block]
            bot.append(mean(bot_33))
            med.append(mean(med_33))
            top.append(mean(top_33))
            all.append(mean(distance_array))
            # fig, axs = plt.subplots(2,2)
            # for row,col,data,title in zip([0,0,1,1],[0,1,0,1],
            #                               [bot_33,med_33,top_33,distance_array],['bottom 33%','med 33%','top 33%','all']):
            #     axs[row,col].hist(data,bins=25,alpha=0.8,color='black',density=False)
            #     axs[row,col].set_title(title)
            #     axs[row,col].set_xlabel('Distance')
            #     axs[row,col].set_ylabel('frequency')
            # fig.tight_layout()
            # plt.savefig('D:/data/native_protein_digestion/10282021/distance_distribution_100prots/%s.png' % prot)
            # plt.close(fig)

    fig, axs = plt.subplots(2,2)
    for row,col,data,title in zip([0,0,1,1],[0,1,0,1],
                                  [bot,med,top,all],['bottom 33%','med 33%','top 33%','all']):
        axs[row,col].hist(data,bins=20,alpha=0.8,color='black',density=False)
        axs[row,col].set_title(title)
        axs[row,col].set_xlabel('Distance')
        axs[row,col].set_ylabel('frequency')
    fig.suptitle('Random 100 proteins average residue distance distribution', fontsize=12)
    plt.show()