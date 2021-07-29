# import __main__
# __main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
import re
import pymol
import numpy as np
from collections import defaultdict
import sys
import argparse
import time


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def peptide_counting(peptide_tsv_file):
    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string,regex_pat='\w{1}\[\d+\.?\d+\]'):
    freq_array = np.zeros(len(protein_seq_string))
    PTM_sites_counting = defaultdict(int)
    PTM_loc_list = []
    #print peptide_list

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation
    for pep, new_pep, PTM in zip(peptide_list, new_pep_list,PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:
            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) -1
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:

                    PTM_index = pep.find(ele)
                   #PTM_site = pep[PTM_index] # single amino acid
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos+PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)
    return freq_array, PTM_loc_list, PTM_sites_counting


def modified_peptide_from_psm(psm_path):
    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]',line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list


def show_cov_3d(peptide_list, protein_seq, pdb_file, png_sava_path):
    """
    show 3d coverage map based on peptide list and protein_seq
    :param peptide_list:
    :param protein_seq:
    :param pdb_file:
    :param png_sava_path:
    :return:
    """
    freq_array, ptm_loc_list, PTM_sites_counting = freq_array_and_PTM_index_generator(peptide_list,protein_seq)
    time_start = time.time()
    # print (ptm_loc_list)
    print ('coverage: %f%%' %(np.count_nonzero(freq_array)/len(freq_array)*100))
    print (f'ptm sites counting: {PTM_sites_counting}')
    # open pdb file with pymol
    pdb_name = pdb_file.split('/')[-1].split('.')[0]
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable()
    # print(pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('white')

    # highlight covered region
    max_freq = np.max(freq_array)
    for i,j in enumerate(freq_array):
        if i not in ptm_loc_list:
            if freq_array[i]==0:
                pymol.cmd.color('grey','resi %i'%(i+1))
            elif 1 <= freq_array[i] < 0.2*max_freq:
                pymol.cmd.color('paleyellow','resi %i'%(i+1))
            elif 0.2*max_freq <= freq_array[i] < 0.4*max_freq:
                pymol.cmd.color('tv_yellow', 'resi %i' % (i + 1))
            elif 0.4*max_freq <= freq_array[i] < 0.6*max_freq:
                pymol.cmd.color('yelloworange', 'resi %i' % (i + 1))
            elif 0.6*max_freq <= freq_array[i] < 0.8*max_freq:
                pymol.cmd.color('tv_orange', 'resi %i' % (i + 1))
            else:
                pymol.cmd.color('sand', 'resi %i' % (i + 1))
        else:
            pymol.cmd.color('red', 'resi %i' % (i + 1))

    pymol.cmd.png(png_sava_path)
    print (f'time used for mapping: {time.time()-time_start}')
    # Get out!
    # pymol.cmd.quit()



if __name__=='__main__':


    """
    pdb_name = pdf_file.split('/')[-1].split('.')[0]
    pymol.finish_launching()
    pymol.cmd.load(pdf_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable()
    print(pymol.cmd.get_names())
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)
    pymol.cmd.bg_color('white')
    # pymol.cmd.color('red', 'ss h')  # coloring helix
    # pymol.cmd.color('yellow', 'ss s')  # coloring b-sheet
    # pymol.cmd.select('toBecolored','resi 1-10') # selection
    for i, j in zip([1, 10, 20, 30], [5, 15, 25, 50]):  # coloring residue range
        pymol.cmd.color('yellow', 'resi %i-%i' % (i, j))
    pymol.cmd.png("%s.png" % (pdb_name))  # save as png
    """
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="coverage map on 3D Alphafold predicted protein structure")
    parser.add_argument("--psm_tsv", required=True, help="psm tsv file path")
    parser.add_argument("--fasta", required=True, help="fasta file path")
    parser.add_argument("--protein_id", required=True, help="protein ID to be mapped")
    parser.add_argument("--pdb", required=True, help="pdb file path")
    parser.add_argument("--output_png", help="PNG output file, if not given, save under same directory")
    args = parser.parse_args()

    psm_tsv, fasta_path, pdb_path, protein_id = args.psm_tsv, args.fasta, args.pdb, args.protein_id
    png_output = args.output_png if args.output_png else '%s.png' % protein_id

    protein_dict = fasta_reader(fasta_path)
    psm_list = modified_peptide_from_psm(psm_tsv)
    # ptm_loc_list = freq_array_and_PTM_index_generator(psm_list,protein_dict[protein_id])[-1]

    show_cov_3d(psm_list,protein_dict[protein_id],pdb_path,png_output)