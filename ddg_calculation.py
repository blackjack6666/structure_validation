from pdb_operation import read_pdb_fasta, pdb_file_reader
import pymol
import os
import numpy as np
from collections import defaultdict


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


if __name__ == '__main__':
    from pdb_operation import complex_pdb_reader, read_pdb_fasta, pdb_cleaner, pdb_file_reader
    from params import aa_dict
    import wget
    import pandas as pd
    from glob import glob
    import pickle as ppp

    # pdb_file = 'C:/tools/Rosetta/rosetta_src_2021.16.61629_bundle/main/source/bin/test/1dkq.pdb'
    # pdb_fasta = 'C:/tools/Rosetta/rosetta_src_2021.16.61629_bundle/main/source/bin/test/rcsb_pdb_1DKQ.fasta'
    # filter_resi_list = find_core_resi(pdb_file, pdb_fasta)
    # mut_file_gen(filter_resi_list, read_pdb_fasta(pdb_fasta),
    #              'C:/tools/Rosetta/rosetta_src_2021.16.61629_bundle/main/source/bin/test/mut_1dkq_081222.txt')

    # print (resi_sasa_dict)
    pdb_test = 'D:/data/pdb/pdb_human_file/5xtc.pdb'
    pdb_clean = 'D:/data/pdb/pdb_human_file/2g9h_chainA_clean.pdb'
    # res_dict, seq = complex_pdb_reader(pdb_test, chain='B')
    # print (res_dict)
    # pdb_cleaner(pdb_test,pdb_clean,chain='A')

    # clean_res_pos_dict = pdb_file_reader(pdb_clean)
    # print(res_dict)
    # print(clean_res_pos_dict)
    # print (seq)

    # download full covered pdbs from PDB.ORG
    df = pd.read_csv('C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv')
    download_dir = 'F:/full_cover_pdbs/'
    # for each in df.pdbchainID.tolist():
    #     pdb = each.split('>')[1].split('_')[0]+'.pdb'
    #     wget.download(url='http://www.pdb.org/pdb/files/'+pdb,out=download_dir)

    ### clean all downloaded pdbs based on chain name from C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv
    for each in df.pdbchainID.tolist():
        pdb = each.split('>')[1].split('_')[0] + '.pdb'
        print(pdb)
        chain = each.split('_')[1]
        pdb_file = download_dir + pdb
        clean_pdb = download_dir + pdb.split('.pdb')[0] + '_' + chain + '_clean.pdb'
        try:
            pdb_cleaner(pdb_file, clean_pdb, chain)
        except:
            continue

    # pdb_cleaner(download_dir+'5xtd.pdb',download_dir+'5xtd_d_clean.pdb','D')

    ## get protein seq dictionary from pdb files
    # pdb_seq_dict = {}
    # for each in df.pdbchainID.tolist():
    #     print (each)
    #     pdb_name = each.split('>')[1].split('_')[0]
    #     pdb = pdb_name+'.pdb'
    #     chain = each.split('_')[1]
    #     pdb_file = download_dir+pdb
    #     pdb_seq = complex_pdb_reader(pdb_file,chain=chain)[1]
    #     pdb_seq_dict[pdb_name+'_'+chain] = pdb_seq
    # print (pdb_seq_dict)
    # ppp.dump(pdb_seq_dict,open(download_dir+'pdb_seq_dict.p','wb'))
