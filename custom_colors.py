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
    from commons import get_unique_peptide
    from background_changed import show_multiple_color
    pdb_path = 'D:/data/native_protein_digestion/pdb_files/pdb6gmh_O00267.ent'

    peptide_path = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/01h_h2o/peptide.tsv'
    peptide_list = peptide_counting(peptide_path)
    protein_dict = fasta_reader('D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta')
    protein = 'Q96I99'
    protein_seq = protein_dict[protein]
    alpha_fold_pdb = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/AF-'+protein+'-F1-model_v1.pdb'

    color_list = [[1.0,0.25,0.2],[0.21,1.0,0.2],[0.94,0.2,1.0],[0.2,0.72,1]]
    time_points = ['01h_h2o','02h_h2o','04h_h2o','20h_h2o']
    # for time,color in zip(time_points,color_list):
    #     peptide_list = peptide_counting('D:/data/native_protein_digestion/10282021/search_result_4miss/h20/'+time+'/peptide.tsv')
    #     png_path = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/' + prote in +'_'+time+'_alphafold.png'
    #     show_cov_3d_custom_c(peptide_list,protein_seq,alpha_fold_pdb,color,png_path)

    ### map unique peptide only to 3d
    unique_pep_dict = get_unique_peptide(['D:/data/native_protein_digestion/10282021/search_result_4miss/h20/'+time+'/peptide.tsv' for time in time_points])
    # for time,color in zip(time_points,color_list):
    #     peptide_list = unique_pep_dict[time]
    #     png_path = 'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/' + protein + '_' + time + '_alphafold_unique.png'
    #     show_cov_3d_custom_c(peptide_list,protein_seq,alpha_fold_pdb,color,png_path)
    show_multiple_color([v for v in unique_pep_dict.values()],protein_seq,alpha_fold_pdb,color_list,
                        'D:/data/native_protein_digestion/10282021/search_result_4miss/h20/different_color_map/'+protein+'_combined.png')

