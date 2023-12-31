"""functions regarding visualizing user-uploaded custom pdb"""

from collections import defaultdict
from params import *
import re
from Bio.PDB import PDBParser
from Bio import SeqIO
from pymol_test import freq_ptm_index_gen_batch_v2
import pymol
import os
import time
from pymol import cmd
from pymol2glmol import parseDistObj, parseObjMol


def seq_reader(pdb_file):
    # reads a pdb file and output residue sequence
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)

    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    seq_dict = {}
    for model in structure:
        for chain in model:
            seq = ''
            for residue in chain:
                if residue.resname in aa_dict:  # discard h2o, small molecules, etc.
                    seq += aa_dict[residue.resname]
            seq_dict[chain.get_id()] = seq

    # one line
    # seq = ''.join([aa_dict[residue.resname] for model in structure for chain in model for residue in chain if
    #                residue.resname in aa_dict])
    return seq_dict


def seq_reader2(pdb_file):
    # use 'SEQRES' header from pdb to get sequence
    seq = ''
    for record in SeqIO.parse(pdb_file, "pdb-seqres"):
        # add sequence from each chain
        seq += str(record.seq)

    # TODO: if there is no SEQRES header in pdb file, "pdb-atom" instead of "pdb-seqres"
    #  can be used in SeqIO.parse to get sequence from atom lines

    return seq


def pdb_resi_atom_mapper(pdb_file):
    """
    for a custom pdb file (non-AlphaFold), map residue index to atom number
    :param pdb_file: .pdb file
    :return: {'residue_number':(atom_start,atom_end), (atom_start, atom_end)}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file)

    # iterate each model, chain, and residue
    chain_resi_atom_dict = {}
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            resi_atom_dict = {}
            res_num = 1  # count residue locaiton from 1
            for residue in chain:
                res_name = residue.get_resname()
                if res_name in aa_dict:  # only take amino acid into account
                    # res_num = residue.get_id()[1]  # residue number
                    atom_num_list = [atom.get_serial_number() for atom in
                                     residue.get_atoms()]  # get all atoms num for residue
                    # take the atom number range
                    resi_atom_dict[res_num] = [atom_num_list[0], atom_num_list[-1]]
                    res_num += 1
                    # res_atom_number_dict[res_num].add('-'.join([chain_id, str(atom_num_list[0]), str(atom_num_list[-1])]))
            chain_resi_atom_dict[chain_id] = resi_atom_dict
    # reorder dictionary
    # clean_dict = {}
    # for res in res_atom_number_dict:
    #     chain_atom_dict = {}
    #     for each in res_atom_number_dict[res]:
    #         chain_atom_num = each.split('-')
    #         chain_atom_dict[chain_atom_num[0]] = '-'.join(chain_atom_num[1:])
    #     clean_dict[res] = chain_atom_dict
    return chain_resi_atom_dict


def show_3d_custom_pdb(pdb_file,
                       psm_list,
                       protein_dict,
                       regex_color_dict=None,
                       png_sava_path=None,
                       base_path=None):
    """

    :param protein_id:
    :param pdb_file:
    :param id_freq_array_dict: returned by freq_ptm_index_gen_batch_v2
    :param id_ptm_idx_dict: returned by freq_ptm_index_gen_batch_v2
    :param png_sava_path:
    :param base_path: html output base path
    :param regex_color_dict {regex: RGB_list}
    :return:
    """
    time_start = time.time()
    # frequency_array = id_freq_array_dict[protein_id]
    # if id_ptm_idx_dict != {}:
    #     ptm_nonzero_idx_dict = id_ptm_idx_dict[protein_id]
    #
    # else:
    #     ptm_nonzero_idx_dict = None
    chain_frequency_array_dict, chain_ptm_nonzero_idx_dict = peptide_mapping_chains(psm_list, protein_dict,
                                                                                    regex_color_dict)

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
    pymol.cmd.bg_color('black')

    if png_sava_path:
        pymol.cmd.png(png_sava_path)

    print(f'image saved to {png_sava_path}')

    # pymol2glmol, convert pdb to pse and visualize through html
    dump_rep_custom_pdb(pdb_name, chain_frequency_array_dict, chain_ptm_nonzero_idx_dict, regex_color_dict, base_path,
                        pdb_file)
    # dump_rep(pdb_name,base_path)
    print(f'time used for mapping: {pdb_name, time.time() - time_start}')
    # Get out!
    pymol.cmd.quit()


def dump_rep_custom_pdb(name, chain_frequency_array_dict, chain_ptm_nonzero_idx_dict, regex_color_dict, base_path,
                        pdb_file):
    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'),
                            os.path.join(os.getcwd(), 'js'))
        except OSError:
            pass

    try:
        cmd.set('pse_export_version', 1.74)
    except:
        pass

    names = cmd.get_session()['names']
    cmd.set('pdb_retain_ids', 1)

    ret = ''
    for obj in names:
        if (obj == None):
            continue
        if (obj[2] == 0):  # not visible
            continue
        if (obj[1] == 0 and obj[4] == 1 and obj[0] == name):
            ret += parseObjMol(obj)
            # print (ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

    pdb_str = cmd.get_pdbstr(name)
    ret += '\n' + color_getter_custom_pdb(chain_frequency_array_dict, chain_ptm_nonzero_idx_dict, regex_color_dict,
                                          pdb_file)
    print(ret)
    cmd.turn('z', 180)
    view = cmd.get_view()
    cmd.turn('z', 180)
    cx = -view[12]
    cy = -view[13]
    cz = -view[14]
    cameraZ = - view[11] - 150
    fov = float(cmd.get("field_of_view"))
    fogStart = float(cmd.get("fog_start"))
    slabNear = view[15] + view[11]
    slabFar = view[16] + view[11]
    ret += "\nview:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % \
           (cx, cy, cz, cameraZ, slabNear, slabFar, fogStart, fov)
    for i in range(9):
        ret += ",%.3f" % view[i]

    bgcolor = cmd.get_setting_tuple('bg_rgb')[1]

    if len(bgcolor) == 1:
        bgcolor = cmd.get_color_tuple(bgcolor[0])

    ret += "\nbgcolor:%02x%02x%02x" % (int(255 * float(bgcolor[0])), \
                                       int(255 * float(bgcolor[1])), int(255 * float(bgcolor[2])))
    if 'PYMOL_GIT_MOD' in os.environ:
        template = open(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'imported.html')).read(). \
            replace("###INCLUDE_PDB_FILE_HERE###", pdb_str). \
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)
    else:
        template = open('imported.html').read(). \
            replace("###INCLUDE_PDB_FILE_HERE###", pdb_str). \
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)

    if base_path:
        f = open(base_path + name + '.html', 'w')
        print(f'html file to {base_path + name}.html')
    else:
        f = open(name.split('-')[1] + '.html', 'w')
        print('html file to %s' % name.split('-')[1] + '.html')
    f.write(template)
    f.close()


def color_getter_custom_pdb(chain_freq_array_dict, chain_ptm_idx_dict, regex_color_dict, pdb_file):
    """
    for custom pdb file, based on numpy array get blocks of zeros and non_zeros, generate color string for GLMOL
    :param chain_freq_array_dict: peptide_mapping_chains
    :param chain_ptm_idx_dict: peptide_mapping_chains
    :param pdb_file: user loaded .pdb file
    :return:
    """
    import re
    import numpy as np
    from collections import defaultdict

    # from pdb_str get element pos to amino acid pos

    chain_amino_ele_pos_dict = pdb_resi_atom_mapper(
        pdb_file)  # {chain:{residue number:[atom start number, atom end number]}}

    defalt = 'color:0.500,0.500,0.500:'  # grey
    covered = 'color:1.000,0.000,0.000:'  # red

    # get index of zeros and nonzeros
    non_zero_index_dict = {chain: np.nonzero(chain_freq_array_dict[chain])[0] for chain in chain_freq_array_dict}
    zero_index_dict = {chain: np.nonzero(chain_freq_array_dict[chain] == 0)[0] for chain in chain_freq_array_dict}

    # get index blocks of zeros and nonzeros
    cov_pos_block = {
        chain: np.split(non_zero_index_dict[chain], np.where(np.diff(non_zero_index_dict[chain]) != 1)[0] + 1)
        for chain in non_zero_index_dict}
    non_cov_pos_block = {chain: np.split(zero_index_dict[chain], np.where(np.diff(zero_index_dict[chain]) != 1)[0] + 1)
                         for chain in zero_index_dict}
    # print(cov_pos_block, non_cov_pos_block)

    # string concatenate
    for chain in non_cov_pos_block:
        all_block = non_cov_pos_block[chain]  # list of numpy arrays
        if len(all_block[0]) > 0:
            defalt += ','.join([str(chain_amino_ele_pos_dict[chain][each[0] + 1][0]) + '-' + str(
                chain_amino_ele_pos_dict[chain][each[-1] + 1][-1])
                                for each in all_block])
            defalt += ','
    for chain in cov_pos_block:
        all_block = cov_pos_block[chain]  # list of numpy arrays
        if len(all_block[0]) > 0:
            # print (all_block)
            covered += ','.join([str(chain_amino_ele_pos_dict[chain][each[0] + 1][0]) + '-' + str(
                chain_amino_ele_pos_dict[chain][each[-1] + 1][-1])
                                 for each in all_block])
            covered += ','

    # defalt += ','.join([amino_ele_pos_dict[each[0] + 1][chain].split('-')[0] + '-' +
    #                     amino_ele_pos_dict[each[-1] + 1][chain].split('-')[-1]
    #                     for each in non_cov_pos_block for chain in amino_ele_pos_dict[each[0] + 1]])
    # covered += ','.join([amino_ele_pos_dict[each[0] + 1][chain].split('-')[0] + '-' +
    #                      amino_ele_pos_dict[each[-1] + 1][chain].split('-')[-1]
    #                      for each in cov_pos_block for chain in amino_ele_pos_dict[each[0] + 1]])

    # ptm color string concatenate
    ptm_color = ''
    if chain_ptm_idx_dict:
        for chain in chain_ptm_idx_dict:
            ptm_idx_dict = chain_ptm_idx_dict[chain]
            for ptm in ptm_idx_dict:
                ptm_color += 'color:' + ','.join(['%.3f' % (each / 256) for each in regex_color_dict[ptm]]) + ':'
                ptm_color += ','.join([str(chain_amino_ele_pos_dict[chain][idx + 1][0]) + '-'
                                       + str(chain_amino_ele_pos_dict[chain][idx + 1][-1])
                                       for idx in ptm_idx_dict[ptm]])
                ptm_color += '\n'

    return defalt + '\n' + covered + '\n' + ptm_color.rstrip('\n')


def peptide_mapping_chains(peptide_list, protein_dict, regex_dict=None):
    """
    map peptides to multiple chains of sequence
    :param protein_dict: single protein sequence dictionary, {chain_A: sequence, chain_B: sequence, ...}
    :param peptide_list: peptide list with ptms
    :param regex_dict: ptm regex color dictionary
    :return:
    """
    chain_freq_array_dict, chain_ptm_index_dict = freq_ptm_index_gen_batch_v2(peptide_list, protein_dict, regex_dict)[
                                                  :2]

    if regex_dict == None:
        return chain_freq_array_dict, None
    else:
        return chain_freq_array_dict, chain_ptm_index_dict


if __name__ == '__main__':
    import numpy as np

    base_path = 'C:/tools/pymol-exporter-0.01/pymol_exporter/'
    pdb_file = r'D:\data\pdb\pdb_human_file/6t3q.pdb'
    chain_seq_dict = seq_reader(pdb_file)
    # prot_seq = seq_reader2(pdb_file)
    print(chain_seq_dict)
    print(len(chain_seq_dict))
    chain_res_atom_dict = pdb_resi_atom_mapper(pdb_file)
    print(chain_res_atom_dict['L'])
    # print (res_atom_dict[1332])
    # print (seq)
    peptide_list = ['GLRPLFE[100]KKSLEDKTE', 'AEIGMSPWQVMLFRKSPQELL[200]CG', 'DFEEI[300]PE']
    regex_dict = {'E\[100\]': [0, 0, 256], 'L\[200\]': [0, 256, 0], 'I\[300\]': [100, 200, 200]}
    # id_freq_array_dict, id_ptm_idx_dict, h = freq_ptm_index_gen_batch_v2(peptide_list, protein_dict, regex_dict)
    show_3d_custom_pdb(pdb_file, peptide_list, chain_seq_dict, regex_dict, base_path=base_path)
