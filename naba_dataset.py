from background_changed import show_multiple_color
from pymol_test import fasta_reader, modified_peptide_from_psm, peptide_counting, freq_array_and_PTM_index_generator, \
    fasta_reader_gene
import numpy as np
import pickle as pk

fasta_file = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_id = 'Q8TER0'
pdb_file = 'D:/data/alphafold_pdb/UP000005640_9606_HUMAN/' + 'AF-' + protein_id + '-F1-model_v1.pdb'
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
psm_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_f_aggregate_peptide_dict_1219.p', 'rb'))
time_points = ['144_' + each + '_aggregate' for each in ['15', '30', '60', '120', '240']]
color_list = [[156, 0, 252], [236, 142, 56], [155, 255, 119], [255, 255, 0], [51, 255, 255]]
peptide_list = [psm_dict['SNED1'][time] for time in time_points]
show_multiple_color(peptide_list, gene_dict['SNED1'], pdb_file, color_list,
                    png_save_path='D:/thesis/figures/SNED1_5timepoints.png')
