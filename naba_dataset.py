from background_changed import show_multiple_color
from pymol_test import fasta_reader, modified_peptide_from_psm, peptide_counting, freq_array_and_PTM_index_generator, \
    fasta_reader_gene
import numpy as np
import pickle as pk
from pdb_operation import pdb_file_reader
from glob import glob
import pickle as pp

fasta_file = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_id = 'Q9JL15'
pdb_file = 'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/' + 'AF-' + protein_id + '-F1-model_v1.pdb'
pdb_file1 = r'D:\data\pdb\pdb_human_file/6ff7.pdb'
base_path = 'C:/tools/pymol-exporter-0.01/pymol_exporter/'
# protein_dict = fasta_reader(fasta_file)
gene_dict = fasta_reader_gene(fasta_file)
"""
time_points = ['30','120','240','1080']
# peptide_list = [peptide_counting('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'D/peptide.tsv')+
#                 peptide_counting('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'F/peptide.tsv') for each in time_points]

psm_list_2d = [modified_peptide_from_psm('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'D/psm.tsv')+
                modified_peptide_from_psm('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'F/psm.tsv') for each in time_points]

# show_multiple_color(peptide_list,protein_dict[protein_id],pdb_file,color_list=
# [[156,0,252],[236,142,56],[155,255,119],[71,182,221]],
#                     png_save_path='D:/data/Naba_deep_matrisome/07232021_secondsearch/3d_structures/'+protein_id+'_SNED_4timepoints.png')

freq_array_2d = [freq_array_and_PTM_index_generator(psm_list, protein_dict[protein_id])[0] for psm_list in psm_list_2d]
freq_sum = np.sum(freq_array_2d,axis=0)
for i,j in zip(time_points,freq_array_2d):
    print (i,'\n',j.tolist())
print ('aggregated','\n', freq_sum.tolist())
"""
psm_dict = pk.load(open('F:/fred_time_lapse/20230508/analysis/gene_f_aggregate_peptide_dict_0509.p', 'rb'))
# print ([each for each in psm_dict['SNED1']])
time_points = ['144_' + each + '_aggregate' for each in ['15', '30', '60', '120', '240']]
# color_list = [[156, 0, 252], [236, 142, 56], [155, 255, 119], [255, 255, 0], [51, 255, 255]]
color_list = [[61, 101, 245], [242, 44, 44], [44, 242, 44], [186, 44, 242], [245, 136, 20]]
# time_points,color_list = ['144_1080'], [[255,0,0]]

peptide_list = [psm_dict['Lgals8'][time] for time in time_points]
show_multiple_color(peptide_list, gene_dict['Lgals8'], pdb_file1, color_list,
                    base_path=base_path)

"""
pdb_path = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/'
pdb_files = glob(pdb_path + '*F1*.pdb')

new_dict = {}
for each_pdb in pdb_files:
    residue_xyz_dict = pdb_file_reader(each_pdb)[0]
    uniprot_id = each_pdb.split('\\')[-1].split('-')[1]
    new_dict[uniprot_id] = residue_xyz_dict
    print (uniprot_id)
pp.dump(new_dict,open('D:/data/alphafold_pdb/uniprot_residue_xyz_dict_of_dict.p','wb'))
"""
