"""
generate a pymol session file
"""
import pymol
import os
import numpy as np
import time


def pse_gen(protein_id,
            pdb_file,
            id_freq_array_dict,
            id_ptm_idx_dict,
            save_path,
            regex_color_dict=None
            ):
    """
    .pse file is openable by Pymol
    :param protein_id:
    :param pdb_file:
    :param id_freq_array_dict: returned by freq_ptm_index_gen_batch_v2
    :param id_ptm_idx_dict: returned by freq_ptm_index_gen_batch_v2
    :param save_path: pse file save path
    :param regex_color_dict {regex: RGB_list}
    :return:
    """

    frequency_array = id_freq_array_dict[protein_id]
    if id_ptm_idx_dict != {}:
        ptm_nonzero_idx_dict = id_ptm_idx_dict[protein_id]

    else:
        ptm_nonzero_idx_dict = None

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

    # set customized ptm color
    if regex_color_dict:
        for i in regex_color_dict:
            pymol.cmd.set_color(i, regex_color_dict[i])

    # color mapped residues
    for i, j in enumerate(frequency_array):

        if frequency_array[i] == 0:
            pymol.cmd.color('grey', 'resi %i' % (i + 1))

        else:
            # print(i)
            pymol.cmd.color('red', 'resi %i' % (i + 1))

    # color map PTMs
    if ptm_nonzero_idx_dict:
        for ptm in ptm_nonzero_idx_dict:
            for ptm_idx in ptm_nonzero_idx_dict[ptm]:
                pymol.cmd.color(ptm, 'resi %i' % (ptm_idx + 1))
    pymol.cmd.save(save_path)

    pymol.cmd.quit()
