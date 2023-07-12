from pymol_test import freq_array_and_PTM_index_generator,freq_ptm_index_gen_batch_v2
import pymol
from pymol2glmol import *
import numpy as np
from pymol_test import peptide_counting, fasta_reader, modified_peptide_from_psm
from glob import glob


def show_cov_3d(peptide_list,
                protein_seq,
                pdb_file,
                png_sava_path=None,
                base_path=None,
                plot_KR=False,
                color=False,
                alpha=False):
    """
    show 3d coverage map based on peptide list and a protein_seq
    output a png file
    :param peptide_list:
    :param protein_seq:
    :param pdb_file:
    :param png_sava_path:
    :return:
    """

    import time
    time_start = time.time()
    if not plot_KR:
        freq_array, ptm_loc_list, PTM_sites_counting = freq_array_and_PTM_index_generator(peptide_list, protein_seq)
    else:
        import re
        freq_array = np.zeros(len(protein_seq))
        kr_index = [m.end() - 1 for m in re.finditer(r'K(?=[^P])|R(?=[^P])', protein_seq)]
        for each in kr_index:
            freq_array[each] += 1

    # print (freq_array)
    # print (ptm_loc_list)
    # print (f'ptm sites counting: {PTM_sites_counting}')
    # open pdb file with pymol
    pdb_name = os.path.split(pdb_file)[1]
    print (pdb_name)
    # pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable()
    # print(pymol.cmd.get_names())
    # pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('white')
    pymol.cmd.remove('solvent') # optional

    # set customized uncovered color
    if color == False:
        color_input = 'grey'
    else:
        pymol.cmd.set_color('cust_c', color)
        color_input = 'cust_c'

    # highlight covered region
    max_freq = np.max(freq_array)

    # when no transparency is needed
    if not alpha:
        for i, j in enumerate(freq_array):

            if freq_array[i] == 0:
                pymol.cmd.color(color_input, 'resi %i' % (i + 1))

            else:
                # print(i)
                pymol.cmd.color('red', 'resi %i' % (i + 1))
    # when need to set alpha transparency
    else:
        for i, j in enumerate(freq_array):

            if freq_array[i] == 0:
                pymol.cmd.color(color_input, 'resi %i' % (i + 1))
                pymol.cmd.set('cartoon_transparency', alpha, 'resi %i' % (i + 1))
            else:
                # print(i)
                pymol.cmd.color('red', 'resi %i' % (i + 1))
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


def show_multiple_color(psm_list_2d,
                        protein_seq,
                        pdb_file,
                        color_list,
                        ptm_color_list=None,
                        png_save_path=None,
                        base_path=None):
    """
    different time points data on a single 3d model
    :param psm_list_2d: 2D peptide list containing peptides identified from different time point
    :param protein_seq: target protein sequence
    :param pdb_file: pdb file absolute path
    :param color_list: 2d_RGB_list
    :param ptm_color_list: 2d_ptm_RGB_list
    :param png_save_path:
    :return:
    """
    import time
    time_start = time.time()

    # load variables
    freq_array_2d = [(freq_array_and_PTM_index_generator(psm_list, protein_seq)[:2]) for psm_list in psm_list_2d]
    freq_array_clean, ptm_loc_2d = zip(*freq_array_2d)
    # ptm_loc_2d = [ptm_loc_2d[0], [each for each in ptm_loc_2d[1] if each not in ptm_loc_2d[0]]]  # special case, delete for other use
    print (freq_array_clean,ptm_loc_2d)
    color_dict = {str(i)+'_color':j for i,j in zip(range(len(freq_array_2d)), color_list)} # color correspond to sample
    if ptm_color_list:
        ptm_color_dict = {str(i) + '_ptm_color': j for i, j in zip(range(len(freq_array_2d)), ptm_color_list)}

    # initialize pdb file in pymol api
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
    pymol.cmd.bg_color('white')

    # set customized color
    for sample in color_dict:

        pymol.cmd.set_color(sample,color_dict[sample])

    if ptm_color_list:
        for ptm in ptm_color_dict:
            pymol.cmd.set_color(ptm, ptm_color_dict[ptm])

    # color mapped amino acids
    for i in range(len(protein_seq)):
        cov = False
        for j in range(len(freq_array_2d)):
            if freq_array_clean[j][i] == 0:
                continue
            else:
                cov = True
                pymol.cmd.color(str(j)+'_color','resi %i' % (i + 1))

                break
        if cov == False:
            pymol.cmd.color('grey', 'resi %i' % (i + 1))
        ### set transparency on specific residue
        # pymol.cmd.set('cartoon_transparency', 0.8, 'resi %i' % (i + 1))
    # color map PTMs
    if ptm_color_list:
        for ind, val in enumerate(ptm_loc_2d):
            for ptm_loc in val:
                pymol.cmd.color(str(ind) + '_ptm_color', 'resi %i' % (ptm_loc + 1))
                pymol.cmd.set('cartoon_transparency', 0, 'resi %i' % (ptm_loc + 1))  # set ptm no transparent
    if png_save_path:
        pymol.cmd.png(png_save_path)

    time.sleep(1)

    print(f'image saved to {png_save_path}')
    if base_path:
        dump_rep(pdb_name, base_path)
    print (f'time used {time.time()-time_start}')
    # pymol.cmd.quit()
    # pymol.cmd.delete(pdb_name)


def show_3d_batch(psm_list, protein_list, protein_dict, pdb_base_path, png_save_path, time_point='1h',
                  glmol_basepath=None):
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

    id_freq_array_dict, id_ptm_idx_dict = freq_ptm_index_gen_batch_v2(psm_list,protein_dict)

    for id in protein_list:
        time_start = time.time()
        pdb_path = pdb_base_path+'AF-'+id+'-F1-model_v1.pdb'
        if os.path.exists(pdb_path): # if local pdb db has this id entry
            pdb_name = os.path.split(pdb_path)[1]
            print(pdb_name)
            pymol.cmd.load(pdb_path, pdb_name)
            pymol.cmd.disable("all")
            pymol.cmd.enable()
            # print(pymol.cmd.get_names())
            pymol.cmd.hide('all')
            pymol.cmd.show('cartoon')
            pymol.cmd.set('ray_opaque_background', 0)
            pymol.cmd.bg_color('white')
            # stick cysteines
            # pymol.cmd.select('cysteines','resn cys')
            # pymol.cmd.show('sticks','cysteines')

            # highlight covered region
            freq_array = id_freq_array_dict[id]
            max_freq = np.max(freq_array)
            for i, j in enumerate(freq_array):

                if freq_array[i] == 0:
                    pymol.cmd.color('grey', 'resi %i' % (i + 1))

                else:
                    pymol.cmd.color('red', 'resi %i' % (i + 1))
            pymol.cmd.png(png_save_path+id+'_'+time_point+'.png')

            # pymol2glmol, convert pdb to pse and visualize through html
            # dump_rep(pdb_name, glmol_basepath)
            print(f'time used for mapping: {pdb_name, time.time() - time_start}')
        else:
            print (f'{id} not exist in pdb database')
            continue
        # pymol.cmd.delete(pdb_name)
        # pymol.cmd.quit()


def show_3d_batch_2(sample_peptide_dict, protein_list, protein_dict, pdb_base_path, png_save_path):
    """
    3d mapping of peptides from multiple samples (time points)
    :param sample_peptide_dict:
    :param protein_list:
    :param protein_dict:
    :param pdb_base_path:
    :param png_save_path:
    :return:
    """
    sample_id_freq_array_dict = {each: freq_ptm_index_gen_batch_v2(sample_peptide_dict[each], protein_dict)[0]
                                 for each in sample_peptide_dict}
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()

    for id in protein_list:
        time_start = time.time()
        pdb_path = pdb_base_path + 'AF-' + id + '-F1-model_v1.pdb'
        if os.path.exists(pdb_path):  # if local pdb db has this id entry
            pdb_name = os.path.split(pdb_path)[1]
            print(pdb_name)
            pymol.cmd.load(pdb_path, pdb_name)
            pymol.cmd.disable("all")
            pymol.cmd.enable()
            # print(pymol.cmd.get_names())
            pymol.cmd.hide('all')
            pymol.cmd.show('cartoon')
            pymol.cmd.set('ray_opaque_background', 0)
            pymol.cmd.bg_color('white')
            for sample in sample_peptide_dict:
                freq_arry = sample_id_freq_array_dict[sample][id]
                # get index of zeros and nonzeros
                non_zero_index = np.nonzero(freq_arry)[0]
                zero_index = np.nonzero(freq_arry == 0)[0]
                # get index blocks of zeros and nonzeros
                cov_pos_block = np.split(non_zero_index, np.where(np.diff(non_zero_index) != 1)[0] + 1)
                non_cov_pos_block = np.split(zero_index, np.where(np.diff(zero_index) != 1)[0] + 1)
                for non_cov in non_cov_pos_block:
                    pymol.cmd.color('grey', 'resi %i-%i' % (non_cov[0] + 1, non_cov[-1] + 1))
                for cov in cov_pos_block:
                    if len(cov) > 0:
                        pymol.cmd.color('red', 'resi %i-%i' % (cov[0] + 1, cov[-1] + 1))
                pymol.cmd.png(png_save_path + id + '_' + sample + '.png')
                time.sleep(0.1)
                print(f'{sample} done mapping...')
                # reset all color to green
                pymol.cmd.color('green', 'resi 1-%i' % (len(freq_arry)))

            pymol.cmd.delete(pdb_name)
        else:
            print(f'{id} not exist in pdb database')
            continue
    pymol.cmd.quit()


if __name__ == '__main__':
    import time
    import pandas as pd
    from commons import get_unique_peptide, protein_tsv_reader
    import pickle as pk

    pdb_file_base_path = 'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/'

    # peptide_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_240D/peptide.tsv'
    # time_point_rep = ['1h']
    # psm_tsv_list = ['D:/data/native_protein_digestion/'+ each + '_1_native/psm.tsv' for each in time_point_rep]
    # print(f'{len(psm_tsv_list)} psm files to read...')

    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
    # peptide_list = peptide_counting(peptide_tsv)
    protein_dict = fasta_reader(fasta_file)

    # checked_protein_file = open('D:/data/pdb/pdb_alpha_exactly_same.txt','r').read().split('\n')
    # pdb_to_check = [each.split('\t') for each in checked_protein_file[:-1]]

    """
    psm_list = [psm for file in psm_tsv_list for psm in modified_peptide_from_psm(file)]

    # protein_list = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F_average.xlsx',index_col=0).index.tolist()
    protein_list = ['P61604']
    protein_dict_sub = {prot: protein_dict[prot] for prot in protein_list}
    base_path = 'C:/tools/pymol-exporter-0.01/pymol_exporter/'  # JS folder path includes JS files required in html

    pdb_path_list = [file for each in protein_list for file in glob(pdb_file_base_path + '*' + each + '*.pdb')]
    # show_3d_batch(psm_list,protein_dict_sub,pdb_file_base_path,glmol_basepath=None)

    show_cov_3d(psm_list, protein_dict_sub['P61604'],
                'D:/data/alphafold_pdb/AF-P61604-F1-model_v1.pdb', png_sava_path='D:/data/alphafold_pdb/native_digest_time_laps/HSPE1_1h.png',
                base_path='D:/data/alphafold_pdb/')

    """

    # protein_list = pd.read_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist.xlsx', index_col=0).index
    # protein_to_check = [(i[0],i[1].split('-')[1]) for i in pdb_to_check if i[1].split('-')[1] in protein_list]

    pdb_base_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
    # base_path = 'D:/data/native_protein_digestion/12072021/control/'
    # folders = glob(base_path + '*/')
    # protein_list_chymo = protein_tsv_reader('F:/native_digestion/chymotrypsin_4_16/search/combined_protein.tsv')
    # protein_list_trypsin = protein_tsv_reader('D:/data/native_protein_digestion/12072021/control/combined_protein.tsv',
    #                                           protein_column=3)
    # protein_list = [each for each in protein_list_chymo if each in protein_list_trypsin]
    # print(len(protein_list))
    #
    # time_points = [each.split('\\')[-2] for each in folders]
    # print (time_points)

    ### aggregate psm
    # psm_dict = {val: [psm for file in [base_path + '/' + each + '/psm.tsv' for each in time_points[:idx + 1]]
    #                   for psm in modified_peptide_from_psm(file)]
    #             for idx, val in enumerate(time_points)}
    # psm_dict = {time:modified_peptide_from_psm(base_path+time+'/psm.tsv') for time in time_points}
    ### unique psm
    # psm_dict = get_unique_peptide(glob(base_path + '/*/peptide.tsv'))
    psm_dict = pk.load(open('F:/native_digestion/01242023/time_points/f_peptides_dict.p', 'rb'))
    protein_list = pd.read_csv(
        r'F:\native_digestion\01242023\analysis/distance_to_center_times_normIntensity_filter_norm01_removelast2.tsv',
        sep='\t', index_col=0).index.tolist()
    # control_df = pd.read_excel('D:/data/native_protein_digestion/12072021/control/spearman_corr_pval_nofill.xlsx',index_col=0)
    # protein_cand_list = control_df.loc[(control_df['spearman correlation']<0)&(control_df['p value']<0.05)].index.tolist()
    time_points = [f for f in psm_dict]
    show_3d_batch_2(psm_dict, protein_list, protein_dict, pdb_base_path,
                    png_save_path='F:/native_digestion/01242023/3d_mapping/')
    # for each_protein in protein_to_check:
    # for each_protein in ['Q14157']:
    # for each_protein in protein_list:
    #     pdb_file_name = 'AF-' + each_protein + '-F1-model_v1.pdb'
    #     if os.path.exists(pdb_base_path + pdb_file_name):
    #         print(each_protein)
    #         for val in time_points:
    #             print(val)
    #             psm_list = psm_dict[val]
    #             show_cov_3d(psm_list, protein_dict[each_protein], pdb_base_path + pdb_file_name,
    #                         png_sava_path='F:/native_digestion/01242023/3d_mapping/' + each_protein + '_' + val + '_trypsin.png')
    #
    #     else:
    #         print(f"{pdb_file_name} not existed")

    ### map peptides to pdbs
    # for each_protein in protein_to_check:
    # for each_protein in ['P30086']:
    #     pdb_file_name = 'D:/data/pdb/robetta_models_PEBP1_P30086.pdb'
    #
    #
    #     for val in time_points:
    #         print(val)
    #         psm_list = psm_dict[val]
    #         show_cov_3d(psm_list, protein_dict[each_protein], pdb_file_name,
    #                     png_sava_path='D:/data/native_protein_digestion/10282021/search_result_4miss/h20/rossetta_mapping/'
    #                                   + each_protein +'_rossetta_'+ val + '.png')

    ### single protein show
    # show_cov_3d(psm_dict['02h_h2o'], protein_dict['Q96AE4'],
    #             pdb_base_path + 'AF-' + 'Q96AE4' + '-F1-model_v1.pdb', alpha=0.6,
    #             png_sava_path='D:/data/SCV/FUBP1_02h_alphafold.png')

    # show_cov_3d(['QQAAYYAQTSPQGMPQHPPAPQGQ'], protein_dict['Q96AE4'],
    #             pdb_base_path + 'AF-' + 'Q96AE4' + '-F1-model_v1.pdb',
    #             png_sava_path='D:/data/chicago_mass_spec_day/FUBP_18h.png', plot_KR=False)

    # for i in ['1h','2h','4h','18h']:
    #     show_cov_3d(psm_dict[i],protein_dict['P61604'],pdb_base_path+'AF-P61604-F1-model_v1.pdb',
    #             'D:/data/native_protein_digestion/dialysis_cassette_cov_png/P61604'+'_'+i+'.png')

    # sub_protein_dict = {each:protein_dict[each] for each in protein_list}
    # show_3d_batch(psm_dict['1h'],sub_protein_dict,pdb_base_path,
    #               png_save_path='D:/data/native_protein_digestion/dialysis_cassette_cov_png/',time_point='1h')
