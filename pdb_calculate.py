"""
calculate some index (distance, density, solvent accessibility area) for uniprot-mapped pdb files
calculate shape of protein structure based on pdb file
"""
from pdb_operation import read_pdb_fasta, pdb_file_reader
import pymol
import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import time


def find_core_resi(pdb_file, pdb_fasta_file):
    """
    find residues inside protein core based on sasa (surface accessibility area)
    :param pdb_file:
    :return: a list with core residue positions (sasa=0)
    """

    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()

    pymol.cmd.set('dot_solvent', 1)
    pymol.cmd.set('dot_density', 3)  # surface area
    pymol.cmd.set('solvent_radius', 3)
    pdb_name = os.path.split(pdb_file)[1]
    pymol.cmd.load(pdb_file, pdb_name)

    protein_seq = read_pdb_fasta(pdb_fasta_file)
    residue_sasa_dict = {}

    for i in range(len(protein_seq)):
        print(i + 1, protein_seq[i])
        residue_sasa_dict[i + 1] = pymol.cmd.get_area(f'resi {i + 1}')
    pymol.cmd.delete(pdb_name)
    print(pdb_file + ' done')
    filter_resi = [each for each in residue_sasa_dict if residue_sasa_dict[each] == 0]

    return filter_resi


def mut_file_gen(resi_pos_list, protein_seq, output_mut: str):
    """
    generate a mutant file needed for Rosetta monomer ddg function
    :param resi_pos_list: residue positions to be mutant
    :param protein_seq: protein aa seq
    :return:
    """
    import random

    hydropobic_aa = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W']
    hydrophilic_aa = ['D', 'E', 'R', 'K', 'H', 'N', 'Q', 'S', 'T', 'Y']
    mutant_dict = defaultdict(set)
    for i in resi_pos_list:
        aa = protein_seq[i - 1]
        if aa in hydrophilic_aa:
            for each in hydropobic_aa:
                mutant_dict[str(i) + '_' + aa].add(each)
            # mutant_dict[str(i)+'_'+aa].add(random.choice(hydropobic_aa))
        else:
            continue

    total_mutant = sum([len(mutant_dict[each]) for each in mutant_dict])
    with open(output_mut, 'w', newline='\n') as f_w:
        f_w.write('total ' + str(total_mutant) + ' \n')
        for each in mutant_dict:
            pos, aa = each.split('_')
            for each_mut in mutant_dict[each]:
                f_w.write('1 \n')
                f_w.write(aa + ' ' + pos + ' ' + each_mut + ' \n')

    return mutant_dict


def structure_volume(pdb_file):
    """
    calculate the structure volume given a pdb file, by projecting onto a 3D grid, and calculate how many grid points
    being occupied
    :param pdb_file:
    :return:
    """
    import itertools
    protein_cood_dict = pdb_file_reader(pdb_file)[0]
    normalized_protein_coord_dict = {}

    # numpy array the atom coordinates
    np_coord_array = np.array([each for each_2d in protein_cood_dict.values() for each in each_2d])

    # get x, y, z range
    x_diff, y_diff, z_diff = [(i - j, j) for i, j in
                              zip(np.max(np_coord_array, axis=0), np.min(np_coord_array, axis=0))]
    # print(x_diff, y_diff, z_diff)
    # coordinate range
    # coord_range = [int(x_diff[0]*10), int(y_diff[0]*10), int(z_diff[0]*10)]
    coord_range = [int(x_diff[0]), int(y_diff[0]), int(z_diff[0])]  # 3d range

    # a 3d grid to project protein structure
    prot_3d_array = np.zeros((coord_range[0] + 1, coord_range[1] + 1, coord_range[2] + 1), dtype=np.int8)

    total_volume = coord_range[0] * coord_range[1] * coord_range[2]

    # normalize x y z coordinates
    for _ in np_coord_array:
        x, y, z = _
        # print(x,y,z)
        new_x, new_y, new_z = int((x - x_diff[1])), int((y - y_diff[1])), int((z - z_diff[1]))
        # print(new_x, new_y, new_z)
        # normalized_protein_coord_dict[_] = (new_x, new_y, new_z)
        prot_3d_array[new_x - 1:new_x + 1, new_y - 1:new_y + 1,
        new_z - 1:new_z + 1] += 1  # also add surrounding of each atom

    volume = np.count_nonzero(prot_3d_array)  # number of grid points occupied by atoms
    density = volume / total_volume  # percentage of protein occupied grid points to total number of grid points in a cube
    # print(volume, density)
    return volume, density


def shape_difference(pdb_file):
    """
    tell the shape of protein structure by pdb coordinates, if it's linearized
    :param pdb_file:
    :return:
    """
    protein_cood_dict = pdb_file_reader(pdb_file)[0]
    normalized_protein_coord_dict = {}

    # numpy array the atom coordinates
    np_coord_array = np.array([each for each_2d in protein_cood_dict.values() for each in each_2d])

    # get x, y, z range
    x_diff, y_diff, z_diff = [(i - j, j) for i, j in
                              zip(np.max(np_coord_array, axis=0), np.min(np_coord_array, axis=0))]
    # print(x_diff, y_diff, z_diff)
    # check the ratio of each axis
    range_divide = [x_diff[0] / y_diff[0], x_diff[0] / z_diff[0], y_diff[0] / z_diff[0]]
    linear = False
    for each in range_divide:
        if each >= 3 or each <= 1 / 3:  # if one axis is 3 times longer than the other, it looks like linear
            linear = True
    return linear


def angle(directions):
    """Return the angle between vectors
    """
    vec2 = directions[1:]
    vec1 = directions[:-1]

    norm1 = np.sqrt((vec1 ** 2).sum(axis=1))
    norm2 = np.sqrt((vec2 ** 2).sum(axis=1))
    cos = (vec1 * vec2).sum(axis=1) / (norm1 * norm2)
    return np.arccos(cos)


def fit_into_function(array_2d, proteinid='XXX', save_path=None):
    """
    fit a distance/sasa array into function and decide number of turning points, first time of turning points
    :param array_2d: array-like trend, time series
    :return:
    """
    from rdp import rdp
    x, y = zip(*array_2d)

    simplified_trajectory = rdp(array_2d, epsilon=1)
    sx, sy = np.array(simplified_trajectory).T
    # Visualize trajectory and its simplified version.
    # Define a minimum angle to treat change in direction
    # as significant (valuable turning point).
    min_angle = np.pi / 3

    # Compute the direction vectors on the simplified_trajectory.
    directions = np.diff(simplified_trajectory, axis=0)
    # print (simplified_trajectory)
    theta = angle(directions)
    # print (theta)
    # Select the index of the points with the greatest theta.
    # Large theta is associated with greatest change in direction.
    idx = np.where(theta > min_angle)[0] + 1
    turning_points = sx[idx], sy[idx]

    # get first turning point (potentially unfolding time point)
    first_turning_point = sx[idx[0]] if len(idx) != 0 else 'monotonic'

    # Visualize valuable turning points on the simplified trjectory.
    if save_path:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y, 'kx-', label='original data')
        ax.plot(sx, sy, 'gx-', label='simplified trajectory')
        ax.plot(sx[idx], sy[idx], 'ro', markersize=7, label='turning points')
        ax.set_xlabel("time point")
        ax.set_ylabel("atom density")
        ax.set_title(proteinid)
        ax.legend(loc='best')
        print(f'saving figure to {save_path + proteinid}atom_dens.png ...')
        plt.savefig(save_path + proteinid + 'atom_dens.png')
        plt.close()
    # plt.show()
    return turning_points, first_turning_point


def batch_analysis_turning_point(df: pd.DataFrame, df_output='turning_point.tsv'):
    """
    correlate first turning point with density of structure
    :param df: pd dataframe with uniprot id as index, time points as columns, values are distance/sasa, etc
    :return:
    """
    np_array = df.to_numpy()
    ind_list = df.index.tolist()
    time_points = np.array([1, 2, 3, 4, 5, 6, 7])
    new_df = pd.DataFrame(index=ind_list, columns=['first turning point'])

    for prot, trend in zip(ind_list, np_array):
        nan_boolearn = np.isnan(trend)
        num_nan = np.count_nonzero(nan_boolearn)  # count number of nan
        if num_nan <= 3:  # if nan data is less than 3
            print(prot)
            non_nan = trend[~nan_boolearn]
            non_nan_idx = np.argwhere(~nan_boolearn).flatten()
            non_nan_times = time_points[non_nan_idx]
            array_2d = list(zip(non_nan_times, non_nan))
            first_turn = fit_into_function(array_2d, prot)[1]
            new_df.loc[prot, 'first turning point'] = first_turn
            # print (num_nan,trend, trend[~nan_boolearn], np.argwhere(~nan_boolearn).flatten())
        else:
            continue
    new_df.to_csv(df_output, sep='\t')


def mobidb_data_clean(data_frame: pd.DataFrame):
    """
    clean data curated from mobidb
    :param data_frame: tsv file from https://mobidb.org/
    :return:
    """
    prot_list = data_frame['acc'].unique()
    prot_disorder_dict = {}
    for prot in prot_list:
        sub_df = data_frame[data_frame['acc'] == prot]

        if 'curated-disorder-priority' in sub_df['feature'].tolist():
            disorder_fraction = float(sub_df[sub_df['feature'] == 'curated-disorder-priority']['content_fraction'])
        elif 'prediction-disorder-priority' in sub_df['feature'].tolist():
            disorder_fraction = float(sub_df[sub_df['feature'] == 'prediction-disorder-priority']['content_fraction'])
        else:
            disorder_fraction = 0
        print(prot, disorder_fraction)
        prot_disorder_dict[prot] = disorder_fraction
    return prot_disorder_dict


if __name__ == '__main__':
    from pdb_operation import complex_pdb_reader, read_pdb_fasta, pdb_cleaner, pdb_file_reader
    from params import aa_dict
    import wget
    import pandas as pd
    from glob import glob
    import pickle as ppp
    import json

    # structure_volume(r'D:\data\alphafold_pdb\UP000005640_9606_HUMAN/AF-P04075-F1-model_v1.pdb')

    # pdb_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
    # protein_list = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx',index_col=0).index
    # pdb_files = [pdb_path + 'AF-' + each + '-F1-model_v1.pdb' for each in protein_list if
    #              os.path.exists(pdb_path + 'AF-' + each + '-F1- .pdb')]
    #
    # for each in pdb_files:
    #
    #     linear = shape_difference(each)
    #     if linear:
    #         print(each.split('-')[1])
    # array2d = [[1,168.5],[2,266.6666667],[3,210.8888889],[4,221.5],[5,np.nan],[6,255.5]]
    # fit_into_function(array2d)

    # density_df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_KR_density_15A.xlsx',index_col=0)
    # batch_analysis_turning_point(density_df,df_output='D:/data/native_protein_digestion/12072021/turning_points/first_turn.tsv')
    # turning_df = pd.read_csv('D:/data/native_protein_digestion/12072021/turning_points/first_turn.tsv',index_col=0,delimiter='\t')
    # turning_df1 = turning_df.dropna()
    # valid_prot_list = turning_df1.index.tolist()
    # url_str = 'https://mobidb.org/api/download?format=tsv&projection=&acc='+'%2C'.join(valid_prot_list)+'&proteome=UP000005640'
    # print (url_str)
    # prot_disorder_dict = mobidb_data_clean(pd.read_csv('D:/data/native_protein_digestion/12072021/control/mobidb_result_2022-11-22T17_27_51.911Z.tsv',sep='\t'))

    # for prot,turn in zip(turning_df1.index, turning_df1['first turning point'].tolist()):
    #     density = structure_volume('D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-'+prot+'-F1-model_v1.pdb')[1]
    # turning_df.loc[prot,'density'] = density
    # disorder_fraction = prot_disorder_dict[prot]
    # turning_df1.loc[prot, 'disorder fraction'] = disorder_fraction

    # turning_df1.to_csv('D:/data/native_protein_digestion/12072021/turning_points/first_turn_disorder_nadrop.tsv',sep='\t')
    disorder_df = pd.read_csv('D:/data/native_protein_digestion/12072021/turning_points/first_turn_disorder_nadrop.tsv',
                              index_col=0, delimiter='\t')
    disorder_df = disorder_df.boxplot(by='first turning point', grid=False, rot=45)
    plt.show()

    # pdb_file = 'C:/tools/Rosetta/rosetta_src_2021.16.61629_bundle/main/source/bin/test/1dkq.pdb'
    # pdb_fasta = 'C:/tools/Rosetta/rosetta_src_2021.16.61629_bundle/main/source/bin/test/rcsb_pdb_1DKQ.fasta'
    # filter_resi_list = find_core_resi(pdb_file, pdb_fasta)
    # mut_file_gen(filter_resi_list, read_pdb_fasta(pdb_fasta),
    #              'C:/tools/Rosetta/rosetta_src_2021.16.61629_bundle/main/source/bin/test/mut_1dkq_081222.txt')

    # print (resi_sasa_dict)
    pdb_test = 'F:/full_cover_pdbs/6hn3.pdb'
    pdb_clean = 'F:/full_cover_pdbs/6hn3_A_clean.pdb'
    # res_dict, seq = complex_pdb_reader(pdb_test, chain='A')
    # print (res_dict)
    # pdb_cleaner(pdb_test,pdb_clean,chain='1')
    # print (seq)
    # new_seq = pdb_file_reader(pdb_clean)[1]
    # print (new_seq)
    # clean_res_pos_dict = pdb_file_reader(pdb_clean)
    # print(res_dict)
    # print(clean_res_pos_dict)
    # print (seq)

    # download full covered pdbs from PDB.ORG
    # df = pd.read_csv('C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv')
    download_dir = 'F:/full_cover_pdbs/'
    # for each in df.pdbchainID.tolist():
    #     pdb = each.split('>')[1].split('_')[0]+'.pdb'
    #     wget.download(url='http://www.pdb.org/pdb/files/'+pdb,out=download_dir)

    ### clean all downloaded pdbs based on chain name from C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv
    # for each in df.pdbchainID.tolist():
    #     pdb = each.split('>')[1].split('_')[0] + '.pdb'
    #     print(pdb)
    #     chain = each.split('_')[1]
    #     pdb_file = download_dir + pdb
    #     clean_pdb = download_dir + pdb.split('.pdb')[0] + '_' + chain + '_clean.pdb'
    #     try:
    #         pdb_cleaner(pdb_file, clean_pdb, chain)
    #     except:
    #         continue

    # pdb_cleaner(download_dir+'5xtd.pdb',download_dir+'5xtd_d_clean.pdb','D')

    ## get protein seq dictionary from pdb files
    # count = 0
    # pdb_seq_dict = {}
    # for each in df.pdbchainID.tolist():
    #
    #     pdb_name = each.split('>')[1].split('_')[0]
    #     pdb = pdb_name+'.pdb'
    #     chain = each.split('_')[1]
    #     # pdb_file = download_dir+pdb
    #     pdb_file = download_dir+ pdb_name+'_'+chain+'_clean.pdb'
    #     print (pdb_file)
    #     # pdb_seq = complex_pdb_reader(pdb_file,chain=chain)[1]
    #     pdb_seq = pdb_file_reader(pdb_file)[1]
    #     pdb_seq_dict[pdb_name+'_'+chain] = pdb_seq
    #     count+=1
    #     print(count)
    # print (pdb_seq_dict)
    # ppp.dump(pdb_seq_dict,open(download_dir+'pdb_seq_dict.p','wb'))
    # pdb_seq_dict = ppp.load(open('F:/full_cover_pdbs/pdb_seq_dict.p', 'rb'))

    ## computation on pdb file, distance, density and sasa
    """
    from commons import protein_tsv_reader, get_unique_peptide
    from pymol_test import mapping_KR_toarray
    import pickle
    from pdb_operation import cov_KR_density, residue_density_cal2, reorder_pdb, sasa_pdb

    df_prot_pdb = pd.read_csv('C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv')
    uniprot_pdb_dict = {prot.split('>')[1]: '_'.join(pdb.split('>')[1].split('_')[:2])
                        for prot, pdb in zip(df_prot_pdb['queryID'], df_prot_pdb['pdbchainID'])}
    #
    pdb_base = 'F:/full_cover_pdbs/'
    protein_tsv = 'D:/data/native_protein_digestion/12072021/control/combined_protein.tsv'
    protein_list = protein_tsv_reader(protein_tsv, protein_column=3)

    print(len([k for k in pdb_seq_dict]))
    sub_protein_dict = {}
    for prot in protein_list:
        if prot in uniprot_pdb_dict:
            pdb_name = uniprot_pdb_dict[prot]
            sub_protein_dict[prot + '_' + pdb_name] = pdb_seq_dict[pdb_name]


    base_path = 'D:/data/native_protein_digestion/12072021/control/'
    folders = [base_path + folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]
    time_points = [each.split('/')[-1] for each in folders]
    pep_path_list = [each + '/peptide.tsv' for each in folders]
    psm_path_list = [each + '/psm.tsv' for each in folders]
    unique_peptide_dict = get_unique_peptide(pep_path_list)
    df = pd.DataFrame(index=protein_list, columns=time_points)  # some protein entry does not have pdb
    # density_dict = ppp.load(open('F:/full_cover_pdbs/mapped_pdb_15A_KR_density.pkl','rb'))
    sasa_dict = json.load(
        open(r'F:\native_digestion\Accessible region\sasa_area15A_dict_pdb.json', 'r'))  # from bowei's data

    for pep_tsv in pep_path_list:
        print(pep_tsv)
        # peptide_list = peptide_counting(pep_tsv)
        peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]
        freq_array_dict, freq_array_index_dict = mapping_KR_toarray(peptide_list, sub_protein_dict)
        for prot in protein_list:
            print(prot)

            if prot in uniprot_pdb_dict:  # if prot has full pdb coverage
                pdb_name = uniprot_pdb_dict[prot]  # pdbid_chain
                pdb_file_path = pdb_base + pdb_name + '_clean.pdb'
                # print (pdb_file_path)
                pdb_seq = sub_protein_dict[prot+'_'+pdb_name]
                # print (pdb_seq)
                # print (pdb_file_path)
                # print (pdb_seq_dict[uniprot_pdb_dict[prot]])
                # plddt_dict = residue_plddt_retrieve(pdb_file_path)

                freq_array = freq_array_dict[prot + '_' + pdb_name]
                # ave_KR_density = cov_KR_density(freq_array, density_dict[pdb_name+'_clean'])
                ave_sasa = cov_KR_density(freq_array, sasa_dict[pdb_name]) if pdb_name in sasa_dict else np.nan
                df.at[prot + '_' + uniprot_pdb_dict[prot], pep_tsv.split('/')[-2]] = ave_sasa

                # else:
                #     print('%s protein len between pdb and fasta is not same' % prot)
            else:
                df.at[prot, pep_tsv.split('/')[-2]] = np.nan
                print(f'{prot} not mapped to pdb')
                continue
    df.to_excel('D:/data/native_protein_digestion/12072021/control/mappdb_KR_sasa_15A.xlsx')
    """

    ### calculation for pdbs
    """
    import multiprocessing
    from glob import glob
    import time

    input_list_tuples = [(download_dir + each + '_clean.pdb', pdb_seq_dict[each]) for each in pdb_seq_dict]
    print(f'in total {len(input_list_tuples)} clean pdb files')

    start = time.time()
    with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as pool:
        result = pool.map(residue_density_cal2, input_list_tuples, chunksize=50)
        pool.close()
        pool.join()
    file_density_dict = {k: v for d in result for k, v in d.items()}

    pickle.dump(file_density_dict, open(
        'F:/full_cover_pdbs/mapped_pdb_15A_KR_density.pkl', 'wb'))
    print(time.time() - start)

    # residue_density_cal2(('F:/full_cover_pdbs/1KWM_B_clean.pdb',pdb_seq_dict['1KWM_B']))
    """
