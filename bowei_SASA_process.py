"""
process the files Bowei generated for solvent accessibility of alphafold structures
"""
import re
from glob import glob
import json
import pandas as pd
import matplotlib.pyplot as plt
import pymannkendall as mk


def extract_sasa(input_f):
    """
    extract residue number and average sasa, could use pandas dataframe.group(), should be quicker
    :param input_f:
    :return:
    """

    res_sasa = {}
    with open(input_f, 'r') as f:
        f_split = f.read().split('\nATOM')[1:]

        res_start = 1
        sasa_total = 0
        atom_count = 0
        for line in f_split:
            sasa = abs(float(line.split('   \t')[1].split('\t')[0]))  # sasa area
            # sasa = abs(float(line.split('\t')[-1]))  # last column:sasa volume, get absolute number
            residue_pos = int(re.findall('(?<![A-Z])\d+', line)[1])

            if residue_pos == res_start:
                sasa_total += sasa
                atom_count += 1
            else:
                # res_sasa[res_start] = sasa_total/atom_count
                res_sasa[res_start] = sasa_total
                res_start = residue_pos
                sasa_total = sasa
                atom_count = 1
    # print (sasa_total,atom_count)
    # res_sasa[res_start] = sasa_total/atom_count # last residue
    res_sasa[res_start] = sasa_total
    return res_sasa


if __name__ == '__main__':

    # files = glob(r'F:\native_digestion\alphafold_pdbs_sasa\*\R30/*.contrib')
    # file_dict = {}
    # for f in files:
    #     print (f.split('\\')[-3])
    #     res_sasa = extract_sasa(f)
    #     file_dict[f.split('\\')[-3]] = res_sasa
    # json.dump(file_dict, open(r'F:\native_digestion\sasa_area30A_total_dict.json','w'))

    ## plot
    df_1 = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa.xlsx', index_col=0)
    df_2 = pd.read_excel('F:/native_digestion/sasa_area_15Atotal.xlsx', index_col=0)
    df_3 = pd.read_excel('F:/native_digestion/sasa_volume_15Atotal.xlsx', index_col=0)
    df_4 = pd.read_excel('F:/native_digestion/sasa_area_2Atotal.xlsx', index_col=0)
    df_5 = pd.read_excel('F:/native_digestion/sasa_area_5Atotal.xlsx', index_col=0)
    df_6 = pd.read_excel('F:/native_digestion/sasa_area_30Atotal.xlsx', index_col=0)

    df_list = [df_1, df_2, df_3, df_4, df_5, df_6]
    label_list = ['Pymol SASA area 15A', 'Liang lab calculated SASA area 15A', 'Liang lab calculated SASA volume 15A',
                  'Liang lab calculated SASA area 2A', 'Liang lab calculated SASA area 5A',
                  'Liang lab calculated SASA area 30A']

    color_list = ['#FF0000', '#FFFF00', '#00FF00', '#0080FF', '#7F00FF', '#808080']
    x = range(1, len(df_1.columns) + 1)
    fig, axs = plt.subplots(1, 1, figsize=(8, 8))
    for df, label, color in zip(df_list, label_list, color_list):
        df_mean = df.mean().tolist()
        print(label, mk.original_test(df_mean))
        axs.plot(x, df_mean, linestyle='-', color=color, linewidth=4, label=label)

        axs.set_xticks(x)
        axs.set_xticklabels(list(df_1.columns), fontsize=12, ha="center", rotation=45)
    axs.legend(loc='upper right')
    plt.show()
