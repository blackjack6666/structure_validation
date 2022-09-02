from pymol_test import freq_array_and_PTM_index_generator,freq_ptm_index_gen_batch_v2
import pymol
from pymol2glmol import *
import numpy as np
from pymol_test import peptide_counting, fasta_reader, modified_peptide_from_psm
from glob import glob


def show_cov_3d_custom_c(peptide_list, protein_seq, pdb_file, custom_color=None, png_sava_path=None, base_path=None):
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
    # pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('white')

    if custom_color==None:
        custom_color = [1.0,0.0,0.0]
    else:
        custom_color = custom_color
    pymol.cmd.set_color('custom_color', custom_color)

    # highlight covered region
    max_freq = np.max(freq_array)
    for i, j in enumerate(freq_array):

        if freq_array[i] == 0:
            pymol.cmd.color('grey', 'resi %i' % (i + 1))

        else:
            pymol.cmd.color('custom_color', 'resi %i' % (i + 1))
    # stick cysteines


    # pymol.cmd.select('cysteines', 'resn cys')
    # pymol.cmd.show('sticks', 'cysteines')
    # pymol.cmd.set('stick_color', 'blue')
    # pymol.cmd.set('stick_radius', 0.7)

    # rotate camera
    # pymol.cmd.rotate('x',180)

    if png_sava_path:
        pymol.cmd.png(png_sava_path)
    time.sleep(1)


    print (f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    # dump_rep(pdb_name,base_path)
    pymol.cmd.delete(pdb_name)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # pymol.cmd.save('new_file.pse')
    # Get out!
    # pymol.cmd.quit()


if __name__=='__main__':
    from pdb_operation import pdb_file_reader
    from commons import get_unique_peptide, dta_reader,dta_reader2,dta_result_analyze
    from background_changed import show_multiple_color

    pdb_path = 'D:/data/native_protein_digestion/pdb_files/pdb6gmh_O00267.ent'

    peptide_path = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/01h_h2o/peptide.tsv'
    # peptide_list = peptide_counting(peptide_path)
    protein_dict = fasta_reader('D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta')
    protein = 'P04406'
    protein_seq = protein_dict[protein]
    alpha_fold_pdb = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-'+protein+'-F1-model_v1.pdb'

    color_list = [[1.0,0.25,0.2],[0.21,1.0,0.2],[0.94,0.2,1.0],[0.2,0.72,1]]
    time_points = ['01h_h2o','02h_h2o','04h_h2o','20h_h2o']
    # for time,color in zip(time_points,color_list):
    #     peptide_list = peptide_counting('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/'+time+'/peptide.tsv')
    #     png_path = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/' + prote in +'_'+time+'_alphafold.png'
    #     show_cov_3d_custom_c(peptide_list,protein_seq,alpha_fold_pdb,color,png_path)
    # peptide_list = peptide_counting('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_1080D/peptide.tsv')+\
    #                peptide_counting('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_1080F/peptide.tsv')
    # show_cov_3d_custom_c(peptide_list,protein_seq,alpha_fold_pdb,[0.01,0.71,0.99], 'D:/data/Naba_deep_matrisome/07232021_secondsearch/3d_structures/SNED1_18h.png')
    ### map unique peptide only to 3d
    from glob import glob
    dta_file = 'D:/data/native_protein_digestion/dimethylation/DTASELECT/030515_control_HL_PBS_chy_corr_2015_04_14_14_31470_DTASelect-filter.txt'
    dta_file_list = glob('D:/data/native_protein_digestion/dimethylation/DTASELECT/GAPDH_heated_then_isotope_labeled_sample_2017_07_25_12_229842_DTASelect-filter.txt')
    file_peptide_dict = dta_reader(dta_file_list)
    peptide_2d_list = [[each[0] for each in file_peptide_dict['GAPDH']],
                       [each[0] for each in file_peptide_dict['HGAPDH']]]
    unique_2d = [[each.replace('K','K[28]') for each in peptide_2d_list[0]],
                 [each.replace('K','K[36]') for each in peptide_2d_list[1]]]
    #
    for each in unique_2d:
        print('\n'.join(each))
    c_list = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    ptm_c_list = [[0.020, 0.980, 0.050], [0.035, 0.145, 0.859]]
    show_multiple_color(unique_2d, protein_seq, alpha_fold_pdb, c_list, ptm_color_list=ptm_c_list)
    # dta_file_list = glob('D:/data/native_protein_digestion/dimethylation/DTASELECT/1_*DTASelect-filter.txt') + \
    #             glob('D:/data/native_protein_digestion/dimethylation/DTASELECT/2_*DTASelect-filter.txt')
    """
    peptide_mass_dict = dta_reader2([dta_file])
    prot_iso_label_dict = dta_result_analyze(peptide_mass_dict)
    for prot in prot_iso_label_dict:

        alpha_pdb_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-'+prot+'-F1-model_v1.pdb'
        if os.path.exists(alpha_pdb_path) and prot in protein_dict:
            prot_seq = protein_dict[prot]
    # unique_pep_dict = get_unique_peptide(['D:/data/native_protein_digestion/10282021/search_result_4miss/h20/'+time+'/peptide.tsv' for time in time_points])
            peptide_2d_list = [list(prot_iso_label_dict[prot][label]) for label in ['light','heavy']]
            unique_pep_2d_list = [[each.replace('K','K[28.00]') for each in peptide_2d_list[0]],
                                  [each.replace('K','K[36.00]') for each in peptide_2d_list[1] if each not in peptide_2d_list[0]]]
    # print (unique_pep_2d_list[1])
    # for time,color in zip(time_points,color_list):
    #     peptide_list = unique_pep_dict[time]
    #     png_path = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/' + protein + '_' + time + '_alphafold_unique.png'
    #     show_cov_3d_custom_c(peptide_list,protein_seq,alpha_fold_pdb,color,png_path)

    # show_multiple_color([v for v in unique_pep_dict.values()],protein_seq,alpha_fold_pdb,color_list,
    #                     'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/'+protein+'_combined.png')
            c_list = [[1.0,0.0,0.0],[1.0,0.0,0.0]]
            ptm_c_list = [[0.835,0.949,0.941],[0.153,0.541,0.506]]
            gapdh_pdb = 'D:/data/native_protein_digestion/dimethylation/1znq.pdb'
            show_multiple_color(unique_pep_2d_list,prot_seq,alpha_pdb_path,c_list,ptm_color_list=ptm_c_list,
                                png_save_path='D:/data/native_protein_digestion/dimethylation/pngs/'+prot+'_heavy_light.png')
    """