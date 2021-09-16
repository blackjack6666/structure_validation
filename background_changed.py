from pymol_test import freq_array_and_PTM_index_generator
import pymol
from pymol2glmol import *
import numpy as np
from pymol_test import peptide_counting, fasta_reader, modified_peptide_from_psm
from glob import glob

def show_cov_3d(peptide_list, protein_seq, pdb_file, png_sava_path=None, base_path=None):
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
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('white')

    # highlight covered region
    max_freq = np.max(freq_array)
    for i, j in enumerate(freq_array):

        if freq_array[i] == 0:
            pymol.cmd.color('grey', 'resi %i' % (i + 1))

        else:
            pymol.cmd.color('red', 'resi %i' % (i + 1))
    # stick cysteines


    # pymol.cmd.select('cysteines', 'resn cys')
    # pymol.cmd.show('sticks', 'cysteines')
    # pymol.cmd.set('stick_color', 'blue')
    # pymol.cmd.set('stick_radius', 0.7)


    if png_sava_path:
        pymol.cmd.png(png_sava_path)

    print (f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    dump_rep(pdb_name,base_path)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # pymol.cmd.save('new_file.pse')
    # Get out!
    # pymol.cmd.quit()


if __name__ == '__main__':
    import time
    import pandas as pd


    pdb_file_base_path = 'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/'

    peptide_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_240D/peptide.tsv'
    time_point_rep = ['1h','2h','4h','18h']
    psm_tsv_list = ['D:/data/native_protein_digestion/'+ each + '_1_native/psm.tsv' for each in time_point_rep]
    print(f'{len(psm_tsv_list)} psm files to read...')
    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
    peptide_list = peptide_counting(peptide_tsv)
    protein_dict = fasta_reader(fasta_file)
    psm_list = [psm for file in psm_tsv_list for psm in modified_peptide_from_psm(file)]

    # protein_list = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_ecm_aggregated_D_F_average.xlsx',index_col=0).index.tolist()
    protein_list = ['P08758']
    protein_dict_sub = {prot: protein_dict[prot] for prot in protein_list}
    base_path = 'C:/tools/pymol-exporter-0.01/pymol_exporter/'  # JS folder path includes JS files required in html

    pdb_path_list = [file for each in protein_list for file in glob(pdb_file_base_path + '*' + each + '*.pdb')]
    # show_3d_batch(psm_list,protein_dict_sub,pdb_file_base_path,glmol_basepath=None)

    show_cov_3d(psm_list, protein_dict_sub['P08758'],
                'D:/data/alphafold_pdb/AF-P08758-F1-model_v1.pdb', png_sava_path='D:/data/alphafold_pdb/native_digest_time_laps/ANXA5_18h.png',
                base_path='D:/data/alphafold_pdb/')