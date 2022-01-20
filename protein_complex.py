import json
import os
import pandas as pd
from glob import glob
from pymol_test import modified_peptide_from_psm
from background_changed import show_multiple_color
from pymol_test import fasta_reader, peptide_counting,freq_ptm_index_gen_batch_v2
import pymol
import numpy as np

def show_cov_3d(protein_id,
                pdb_file,
                id_freq_array_dict,
                png_sava_path=None,
                base_path=None):
    """
    show 3d coverage map based on peptide list and a protein_seq
    :param peptide_list:
    :param protein_seq:
    :param pdb_file:
    :param png_sava_path:
    :return:
    """

    import time
    import numpy as np
    time_start = time.time()
    freq_array = id_freq_array_dict[protein_id]
    # print (freq_array)

    # open pdb file with pymol
    pdb_name = os.path.split(pdb_file)[1]
    print (pdb_name)
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
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

    # highlight covered region
    max_freq = np.max(freq_array)
    for i, j in enumerate(freq_array):

        if freq_array[i] == 0:
            pymol.cmd.color('grey', 'resi %i' % (i + 1))

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



if __name__=="__main__":
    ### read protein complex json file into dictionary

    json_path = 'D:/data/coreComplexes.json'
    f = open(json_path)
    complex_list = json.load(f)
    f.close()

    print (len(complex_list))

    ### filter protein complex into human and no disease
    human_complex_list = [each_dict for each_dict in complex_list
                          if each_dict["Organism"] == "Human" and each_dict["Disease comment"] == None]

    human_complex_dict = {
        each["ComplexName"].replace('/', '_'): (each["subunits(UniProt IDs)"], each["subunits(Gene name)"]) for each in
        human_complex_list}
    # for each in human_complex_dict:
    #     print (each, len(human_complex_dict[each][0].split(';')))

    ### read psms from tsv
    pdb_base_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
    base_path = 'D:/data/native_protein_digestion/12072021/control/'
    folders = glob(base_path+'*/')

    time_points = [each.split('\\')[-2] for each in folders]
    print (time_points)

    ### aggregate psm
    """
    psm_dict = {val:[psm for file in [base_path +'/'+ each + '/peptide.tsv' for each in time_points[:idx+1]]
                        for psm in peptide_counting(file)]
                    for idx, val in enumerate(time_points)}
    """
    fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
    protein_dict = fasta_reader(fasta_file)

    ### generate 3d coverage map
    # color_list = [[1.0,0,0],[0,1.0,0.1]]
    """
    proteins_identified = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_dist_unique.xlsx', index_col=0).index
    complex_base_path = 'D:/data/native_protein_digestion/12072021/control_complex_map/'

    pdb_base_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'

    time_point_id_freq_array_dict = {time:freq_ptm_index_gen_batch_v2(psm_dict[time],protein_dict)[0] for time in psm_dict}

    for each_dict in human_complex_list:
        uniprot_ids = each_dict['subunits(UniProt IDs)'].split(";")
        complex_name = each_dict['ComplexName'].replace('/','_')
        genes = each_dict['subunits(Gene name syn)'].split(';')
        protein_id_gene_dict = {prot:gene for prot,gene in zip(uniprot_ids,genes)}
        if any(prot in proteins_identified for prot in uniprot_ids):
            directory = os.path.join(complex_base_path,complex_name)
            # print (directory)
            # print (os.path.join(directory,'asdf'))
            if not os.path.exists(directory):
                os.makedirs(directory)
                for each_prot in uniprot_ids:
                    if each_prot in proteins_identified:
                        print (each_prot)
                        pdb_file_name = 'AF-' + each_prot + '-F1-model_v1.pdb'
                        if os.path.exists(pdb_base_path + pdb_file_name):
                            for val in time_points:
                                show_cov_3d(each_prot,pdb_base_path+pdb_file_name,time_point_id_freq_array_dict[val],
                                            png_sava_path=os.path.join(directory,each_prot+'_'+protein_id_gene_dict[each_prot]+'_'+val+'.png'))
    """

    ### calculate p val for protein complex, fisher exact test
    """
    import os
    from scipy.stats import fisher_exact

    total_id_protein = 1709  # number of prt identifed
    protein_background = 20378  # number of prot in database
    folder_complex = os.listdir("D:/data/native_protein_digestion/12072021/control_complex_map/")
    # print(folder_complex)

    pval_dict = {}
    for complex in folder_complex:
        if complex in human_complex_dict:
            num_id_protein = int(len(glob("D:/data/native_protein_digestion/12072021/control_complex_map/"+complex+'/*.png'))/7)
            total_prot_in_complex = len(human_complex_dict[complex][0].split(';'))  # from CORUM db
            non_id_protein_in_complex = total_prot_in_complex-num_id_protein
            num_other = total_id_protein-num_id_protein
            num_other_nonid = protein_background-total_id_protein-non_id_protein_in_complex

            table = np.array([[num_id_protein,num_other],[non_id_protein_in_complex,num_other_nonid]])
            p_val = fisher_exact(table,alternative='greater')[1]
            pval_dict[complex] = (p_val,total_prot_in_complex, num_id_protein)
    # print(sorted(pval_dict.items(),key=lambda p:p[1][0],reverse=False))
    
    # write to excel
    df = pd.DataFrame()
    for each in pval_dict:
        for i, j in zip(['pval','# proteins in complex CORUM','# proteins identified'],range(3)):
            df.loc[each,i] = pval_dict[each][j]
    df.to_excel('D:/data/native_protein_digestion/12072021/control_complex_map/fisher_exact_pval.xlsx')
    """
