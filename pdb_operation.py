### some functions to operate on PDB file

import re
from collections import defaultdict
from pymol_test import *
import numpy as np
import time


def pdb_file_reader(pdb_file):
    """
    read pdb file and map xyz coordinates of each residue (alphafold pdbs)
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


def pdb_mutiple_reader(pdb_file_list:list):
    """
    read alphafold pdb files corresponding to one protein, because protein too long have multiple alphafold pdb files
    :param pdb_file_list: a list of alphafold pdbs corresponding to one protein entry
    :return:
    """
    residue_atom_xyz = defaultdict(list)

    for pdb in pdb_file_list:
        with open(pdb, 'r') as f_o:
            f_read = f_o.read()
            res_length = int(f_read.split('\nSEQRES')[0].split('    ')[-1].rstrip(' '))
            file_split = f_read.split('\nATOM')[1:]
        for line in file_split:
            residue_atom_xyz[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append([float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)])


def complex_pdb_reader(pdb_file):
    """
    read more complex pdb file, read each residue and 3 coordinates
    :param pdb_file:
    :return:
    """
    from params import aa_dict,aa_reg_str
    residue_atom_xyz = defaultdict(list)

    info_list = []
    with open(pdb_file,'r') as f_o:
        file_read = f_o.read()
        file_split = file_read.split('\nATOM')[1:]
        last_line = file_split[-1]
        if 'HETATM' in last_line:
            last_line = last_line.split('HETATM')[0]
            # for line in file_split[:-1]:
                # residue_atom_xyz[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append(
                #     [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)])
            for line in file_split[:-1]:
                if re.search(aa_reg_str,line):
                    info_list.append((re.search(aa_reg_str,line).group(0),
                                      [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)]))

            info_list.append((re.search(aa_reg_str,last_line).group(0),
                              [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', last_line)]))
            # residue_atom_xyz[int(re.search('\d+(?=\s+[+-]?\d+\.)', last_line).group())].append(
            #     [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', last_line)])
        else:
            # for line in file_split:
                # residue_atom_xyz[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append(
                #     [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)])
            for line in file_split[:-1]:
                if re.search(aa_reg_str,line):
                    info_list.append((re.search(aa_reg_str,line).group(0),
                                     [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)]))

    return info_list


def plddt_retrieve(alphafold_pdb):
    """
    get pLDDT value from alphafold pdb file
    :param alphafold_pdb:
    :return:
    """
    with open(alphafold_pdb,'r') as f_o:
        file_split = f_o.read().split('\nATOM')[1:]

        return [float(re.search('\d+\.\d+(?=\s+[A-Z])',line).group()) for line in file_split]


def residue_plddt_retrieve(alphafold_pdb):
    """
    get pLDDT value for each residue number from alphafold pdb
    :param alphafold_pdb:
    :return:
    """
    with open(alphafold_pdb,'r') as f_o:
        file_split = f_o.read().split('\nATOM')[1:]

    return {int(re.search('\d+(?=\s+[+-]?\d+\.)',line).group()):
                          float(re.search('\d+\.\d+(?=\s+[A-Z])',line).group()) for line in file_split}


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
        # ave_dist = 0
        non_zero_index = np.nonzero(freq_array)[0]
        ave_dist = sum([residue_dist_dict[i+1] for i in non_zero_index])/num_nonzeros
        # for i in non_zero_index:
        #     ave_dist += residue_dist_dict[i+1]

        return ave_dist


def cov_plddt(freq_array,plddt_dict):
    """
    calculate average plddt of covered region in one protein
    :param freq_array:
    :param plddt_dict:
    :return:
    """
    num_nonzeros = np.count_nonzero(freq_array)
    if num_nonzeros == 0:
        return None
    else:
        # ave_dist = 0
        non_zero_index = np.nonzero(freq_array)[0]
        ave_plddt = sum([plddt_dict[i + 1] for i in non_zero_index]) / num_nonzeros
        # for i in non_zero_index:
        #     ave_dist += residue_dist_dict[i+1]

        return ave_plddt


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


def residue_density_cal(input_tuple):
    """
    calculate the number of atoms within a certain range of one residue
    :param alphafold_pdb_file: the alphafold pdb file
    :return:
    """
    time_start = time.time()

    alphafold_pdb_file,protein_seq = input_tuple
    k_density_dict, r_density_dict = {}, {}

    k_index = [m.end() for m in re.finditer(r'K(?=[^P])', protein_seq)]
    r_index = [m.end() for m in re.finditer(r'R(?=[^P])', protein_seq)]

    residue_density_dict = {}
    residue_atom_coord_dict = pdb_file_reader(alphafold_pdb_file)
    # print (residue_atom_coord_dict)
    # residue_xyz_dict = {each:np.mean(residue_atom_coord_dict[each],axis=0) for each in residue_atom_coord_dict}
    xyz_nparray = [v for v in residue_atom_coord_dict.values()]

    # xyz max minus xyz min
    # xyz_range = np.amax(xyz_nparray,axis=0)-np.amin(xyz_nparray,axis=0)

    # check 1/10 of xyz_range for each residue
    # check_range_xyz = xyz_range/10/2

    radius_power2 = 225 # 15A^2 trypsin radius= 1.5 nm
    # for each in residue_atom_coord_dict:
    #
    #     ref = residue_atom_coord_dict[each][-1] # only count C terminal of one residue as ref
    #     # print (ref)
    #
    #     # find CNOS atoms within the radius range
    #     bool_array = [inSphere(j,ref,radius_power2) for i in xyz_nparray for j in i]
    #
    #     # number of CNO atoms within a range for one residue C terminus = number of boolean true - len()
    #     # num_resi_inrange = np.count_nonzero(bool_array)-len(residue_atom_coord_dict[each])
    #     num_resi_inrange = np.count_nonzero(bool_array)
    #     residue_density_dict[each]=num_resi_inrange

    ### get # of density within a range of K/R
    for k in k_index:
        ref = residue_atom_coord_dict[k][-1]
        bool_array = [inSphere(j,ref,radius_power2) for i in xyz_nparray for j in i]
        num_resi_inrange = np.count_nonzero(bool_array)
        k_density_dict[k] = num_resi_inrange
    for r in r_index:
        ref = residue_atom_coord_dict[r][-1]
        bool_array = [inSphere(j,ref,radius_power2) for i in xyz_nparray for j in i]
        num_resi_inrange = np.count_nonzero(bool_array)
        r_density_dict[r] = num_resi_inrange

    # print (time.time()-time_start)
    return {alphafold_pdb_file.split('\\')[-1].split('-')[1]:(k_density_dict,r_density_dict)}


def cov_KR_density(mapped_KR_array,KR_index_density_tuple):
    """
    calculate the average covered K/R density
    :param: mapped_KR_array: mapping of start/end for each peptide on a numpy zero array,same length as protein
    :param: KR_index_density_dict: return by residue_density_cal['proteinid']
    :return:
    """
    k_density_dict,r_density_dict = KR_index_density_tuple
    combined_density_dict = k_density_dict | r_density_dict

    num_nonzeros = np.count_nonzero(mapped_KR_array)
    if num_nonzeros == 0:
        return None
    else:
        non_zero_index = np.nonzero(mapped_KR_array)[0]
        sum_density = 0
        for i in non_zero_index:
            if i+1 in combined_density_dict:
                sum_density+=combined_density_dict[i+1]
            else:
                print(i)

        return sum_density/num_nonzeros


def inSphere(point, ref, radius_power2):
    diff = np.subtract(point, ref)
    # print (diff)
    dist = np.sum(np.power(diff, 2))
    # print (dist)
    # If dist is less than radius^2, return True, else return False
    return dist < radius_power2


def read_pdb_fasta(pdb_fasta):
    """
    read residue sequence from pdb fasta
    :param pdb_fasta:
    :return:
    """
    with open(pdb_fasta, 'r') as f:
        f_split = f.read().split('>')[1:]
        sequence = ''.join([each.split('\n')[-2] for each in f_split])
    return sequence


if __name__ == '__main__':
    import time
    import pickle
    import matplotlib.pyplot as plt
    pdb_file = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-Q9H2X0-F1-model_v1.pdb'
    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
    protein_dict = fasta_reader(fasta_file)
    alphafold_protein_dict = pickle.load(open('D:/data/alphafold_pdb/human_alpha_seq_dict.p','rb'))

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
    ### get unique peptide dict
    pdb_base = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
    """
    from commons import get_unique_peptide

    def protein_tsv_reader(protein_tsv_file):
        with open(protein_tsv_file, 'r') as file_open:
            next(file_open)
            return [line.split("\t")[3] for line in file_open]

    protein_tsv = 'D:/data/native_protein_digestion/12072021/control/combined_protein.tsv'
    protein_list = protein_tsv_reader(protein_tsv)
    sub_protein_dict = {prot:protein_dict[prot] for prot in protein_list}

    
    base_path = 'D:/data/native_protein_digestion/12072021/control/'
    folders = [base_path + folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]
    time_points = [each.split('/')[-1] for each in folders]
    psm_path_list = [each + '/peptide.tsv' for each in folders]
    unique_peptide_dict = get_unique_peptide(psm_path_list)
    print(f'{len(psm_path_list)} psm files to read...')
    """
    ### calculate covered distance/average pLDDT and write to excel
    """
    import pandas as pd
    df = pd.DataFrame(index=protein_list, columns=time_points)  # some protein entry does not have pdb

    for pep_tsv in psm_path_list:
        print (pep_tsv)
        # peptide_list = peptide_counting(pep_tsv)
        peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]
        freq_array_dict = freq_ptm_index_gen_batch_v2(peptide_list,protein_dict)[0]
        for prot in protein_list:
            pdb_file_path = pdb_base+'AF-'+prot+'-F1-model_v1.pdb'
            if os.path.exists(pdb_file_path):
                residue_dist_dict = residue_distance(pdb_file_reader(pdb_file_path))
                # plddt_dict = residue_plddt_retrieve(pdb_file_path)
                if len(residue_dist_dict) == len(protein_dict[prot]):  # filter out those really long proteins
                # if len(plddt_dict) == len(protein_dict[prot]):
                    freq_array = freq_array_dict[prot]
                    cov_dist = cov_distance(freq_array,residue_dist_dict)
                    # ave_cov_plddt = cov_plddt(freq_array,plddt_dict)
                    df.at[prot,pep_tsv.split('/')[-2]] = cov_dist
                    # df.at[prot,pep_tsv.split('/')[-2]] = ave_cov_plddt
                else:
                    print ('%s protein len between pdb and fasta is not same' % prot)
            else:
                continue
    df.to_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx')
    """
    """
    
    from statistics import mean
    import random
    # prots_tocheck = random.sample(protein_list,100)
    prots_tocheck = [each.split('\\')[-1].split('.png')[0] for each in glob('D:/data/native_protein_digestion/10282021/protein_centroid/*')]
    print (prots_tocheck)
    """

    ### calculate covered K/R density and write to excel
    """
    import pandas as pd
    from pymol_test import mapping_KR_toarray

    df = pd.DataFrame(index=protein_list, columns=time_points)  # some protein entry does not have pdb
    KR_density_alpha_dict = pickle.load(open('D:/data/alphafold_pdb/human_file_KR_density_dict.pkl','rb'))
    for pep_tsv in psm_path_list:
        print(pep_tsv)
        # peptide_list = peptide_counting(pep_tsv)
        peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]
        freq_array_dict = mapping_KR_toarray(peptide_list, sub_protein_dict)
        for prot in protein_list:
            print (prot)
            pdb_file_path = pdb_base + 'AF-' + prot + '-F1-model_v1.pdb'
            if os.path.exists(pdb_file_path):
                # residue_dist_dict = residue_distance(pdb_file_reader(pdb_file_path))
                plddt_dict = residue_plddt_retrieve(pdb_file_path)

                # if len(residue_dist_dict) == len(protein_dict[prot]):  # filter out those really long proteins
                if len(plddt_dict) == len(protein_dict[prot]):
                    freq_array = freq_array_dict[prot]
                    ave_KR_density = cov_KR_density(freq_array,KR_density_alpha_dict[prot])
                    # df.at[prot,pep_tsv.split('/')[-2]] = cov_dist
                    df.at[prot, pep_tsv.split('/')[-2]] = ave_KR_density
                else:
                    print('%s protein len between pdb and fasta is not same' % prot)
            else:
                continue
    df.to_excel('D:/data/native_protein_digestion/12072021/control/cov_KR_density.xlsx')
    """

    ### plot 3d and centroid
    from matplotlib import animation


    def rotate(angle):
        ax.view_init(azim=angle)


    # for prot in prots_tocheck:
    for prot in ['Q13148']:
        pdb_file_path = pdb_base+'AF-'+prot+'-F1-model_v1.pdb'

        if os.path.exists(pdb_file_path):
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            print (prot)

            residue_atom_xyz = pdb_file_reader(pdb_file_path)
            # finding coordinates of K and R
            k_ind, r_ind = [m.end() for m in re.finditer(r'K(?=[^P])', protein_dict[prot])], \
                           [m.end() for m in re.finditer(r'R(?=[^P])', protein_dict[prot])]

            k_x, k_y, k_z = zip(*[each for ind in k_ind for each in residue_atom_xyz[ind]])
            r_x, r_y, r_z = zip(*[each for ind in r_ind for each in residue_atom_xyz[ind]])

            # plot all atoms
            # xyz = [each for v in residue_atom_xyz.values() for each in v]

            # get xyz of atoms other than K and R
            xyz = [each for ind in residue_atom_xyz if ind not in k_ind + r_ind for each in residue_atom_xyz[ind]]
            x,y,z = zip(*xyz)

            ax.scatter(x,y,z,marker='o',s=0.5)

            # highlght K and R
            ax.scatter(k_x, k_y, k_z, marker='o', s=0.5, color='green')
            ax.scatter(r_x, r_y, r_z, marker='o', s=0.5, color='orange')

            centroid = find_centroid(residue_atom_xyz)
            # plot centroid
            ax.scatter([centroid[0]],[centroid[1]],[centroid[2]], marker='o', s=8,color='r')
            ax.text2D(0.05, 0.95, "centroid coordinates: %.2f,%.2f,%.2f" % (centroid[0] ,centroid[1],centroid[2]), transform=ax.transAxes)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')

            # Get rid of colored axes planes
            # First remove fill
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            # Now set color to white (or whatever is "invisible")
            # ax.xaxis.pane.set_edgecolor('w')
            # ax.yaxis.pane.set_edgecolor('w')
            # ax.zaxis.pane.set_edgecolor('w')

            # Bonus: To get rid of the grid as well:
            ax.grid(False)
            # make gif
            rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
            rot_animation.save('D:/data/native_protein_digestion/10282021/protein_centroid_median/%s.gif' % prot,
                               dpi=300, writer='imagemagick')
            # plt.show()
            # plt.savefig('D:/data/native_protein_digestion/10282021/protein_centroid_median/%s.png' % prot, dpi=300)


    ### plot residue distance distribution
    """
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
    """

    ### extract pLDDT from all human alphafold pdbs
    import pickle
    from glob import glob
    pdb_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
    pdb_files = glob(pdb_path + '*F1*.pdb')
    # input_list_tuples = [(pdb, alphafold_protein_dict[pdb.split('\\')[-1]]) for pdb in pdb_files]
    # print (len(input_list_tuples))
    """   
    count = 0
    # total_array = []
    pdb_plddt_dict = {}
    for each_pdb in pdb_files:
        count += 1
        print(count)
        pdb_plddt_dict[each_pdb.split('/')[-1].split('-')[1]] = plddt_retrieve(each_pdb)
        # total_array.append(plddt_retrieve(each_pdb))
    # pickle.dump(total_array,open('D:/data/alphafold_pdb/pLDDT_human_2d.pkl','wb'))
    pickle.dump(pdb_plddt_dict,open('D:/data/alphafold_pdb/pLDDT_human_dict.pkl','wb'))
    """

    ### calculate residue density for each alphafold pdb, using multiple cpu cores
    import multiprocessing
    # start = time.time()
    # with multiprocessing.Pool(multiprocessing.cpu_count()-2) as pool:
    #     result = pool.map(residue_density_cal,input_list_tuples,chunksize=500)
    #     pool.close()
    #     pool.join()
    # file_density_dict = {k:v for d in result for k, v in d.items()}
    #
    # pickle.dump(file_density_dict,open('D:/data/alphafold_pdb/human_file_KR_density_dict.pkl','wb'))
    # print (time.time()-start)

    # k_r_density_dict = pickle.load(open('D:/data/alphafold_pdb/human_file_KR_density_dict.pkl','rb'))
    # print (k_r_density_dict['Q8IXR9'])

    # KR_mapped_dict = mapping_KR_toarray(unique_peptide_dict[psm_path_list[1].split('/')[-2]],protein_dict)

