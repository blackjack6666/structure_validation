'''
Pymol to GLmol exporter

Written by biochem_fan, 2011

Modification has been made by Xinhao
'''

from pymol import cmd
from math import cos, sin, pi, sqrt, acos, asin, atan2
import os


def compactSeq(seq):
    seq.sort()
    ret = []
    prev = -9999
    start = -9999
    seq.append(-1)
    i = 0
    while (i < len(seq) - 1):
        if (start >= 0):
            if (seq[i] + 1 != seq[i + 1]):
                ret.append("%d-%d" % (start, seq[i]))
                start = -999
        else:
            if (seq[i] + 1 != seq[i + 1]):
                start = -999
                ret.append(str(seq[i]))
            else:
                start = seq[i]
        i += 1
    return ','.join(ret)


def parseObjMol(obj):
    name = obj[0]
    ids = []
    sphere = []
    trace = []
    ribbon = []
    stick = []
    surface = []
    line = []
    cross = []
    smallSphere = []
    helix = []
    sheet = []
    colors = {}
    for atom in obj[5][7]:
        rep = atom[20] + [0] * 12
        serial = atom[22]
        ss = atom[10]
        bonded = (atom[25] == 1)
        if (rep[5] == 1):
            ribbon.append(serial)
        if (rep[1] == 1):
            sphere.append(serial)
        if (rep[2] == 1):
            surface.append(serial)
        if (rep[7] == 1):
            line.append(serial)
        if (rep[6] == 1):
            trace.append(serial)
        if (rep[4] == 1 and not bonded):
            smallSphere.append(serial)
        if (rep[11] == 1 and not bonded):
            cross.append(serial)
        if (rep[0] == 1 and bonded):
            stick.append(serial)
        if (ss == 'S'):
            sheet.append(serial)
        if (ss == 'H'):
            helix.append(serial)

    # ## retrive color from pymol instance
    #     c = cmd.get_color_tuple(atom[21])
    #     if (not c in colors):
    #         colors[c] = []
    #     colors[c].append(serial)
    #     ids.append("ID %d is %s in resi %s %s at chain %s"\
    #                % (atom[22], atom[6], atom[3], atom[5], atom[1]))
    #
    # for c in colors.keys():  # TODO: better compression
    #     colors[c] = compactSeq(colors[c])

    ret = ''
    ret += "\nsheet:" + compactSeq(sheet)
    ret += "\nhelix:" + compactSeq(helix)
    ret += "\nsurface:" + compactSeq(surface)
    ret += "\nsphere:" + compactSeq(sphere)
    ret += "\ntrace:" + compactSeq(trace)
    ret += "\nribbon:" + compactSeq(ribbon)
    ret += "\nstick:" + compactSeq(stick)
    ret += "\nline:" + compactSeq(line)
    ret += "\nsmallSphere:" + compactSeq(smallSphere)
    ret += "\ncross:" + compactSeq(cross)
    # for c in colors.keys():
    #     ret += "\ncolor:%.3f,%.3f,%.3f:%s" % (c[0], c[1], c[2], colors[c])
    return ret


def parseDistObj(obj):
    if (obj[5][0][3][10] != 1):  # 'show dashed' flag
        return ""
    N = obj[5][2][0][0]
    points = obj[5][2][0][1]
    ret = []
    for p in points:
        ret.append("%.3f" % p)
    color = cmd.get_color_tuple(obj[5][0][2])
    return "\ndists:%.3f,%.3f,%.3f:" % color + ','.join(ret)


def dump_rep(name, base_path=None):
    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'), os.path.join(os.getcwd(), 'js'))
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
            print (ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

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
        template = open(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'imported.html')).read().\
            replace("###INCLUDE_PDB_FILE_HERE###", cmd.get_pdbstr(name)).\
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)
    else:
        template = open('imported.html').read().\
            replace("###INCLUDE_PDB_FILE_HERE###", cmd.get_pdbstr(name)).\
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)

    if base_path:
        f = open(base_path+name + '.html', 'w')
        print (f'html file to {base_path+name}.html')
    else:
        f = open(name.split('-')[1] + '.html', 'w')
        print ('html file to %s' % name.split('-')[1]+'.html')
    f.write(template)
    f.close()


### new functions
def color_getter(freq_array, ptm_idx_dict, regex_color_dict, pdb_str):
    """
    based on numpy array get blocks of zeros and non_zeros, generate color string for GLMOL
    :param freq_array: numpy array
    :param pdb_str: cmd.get_pdbstr()
    :return:
    """
    import re
    import numpy as np
    from collections import defaultdict

    # from pdb_str get element pos to amino acid pos
    pdb_str = pdb_str.split('\nTER')[0].split('\n')
    amino_ele_pos_dict = defaultdict(list)
    for line in pdb_str:
        amino_ele_pos_dict[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append(int(re.search('\d+', line).group()))
    amino_ele_pos_dict = {each:sorted(amino_ele_pos_dict[each]) for each in amino_ele_pos_dict}

    defalt = 'color:0.500,0.500,0.500:'  # grey
    covered = 'color:1.000,0.000,0.000:'

    # get index of zeros and nonzeros
    non_zero_index = np.nonzero(freq_array)[0]
    zero_index = np.nonzero(freq_array==0)[0]

    # get index blocks of zeros and nonzeros
    cov_pos_block = np.split(non_zero_index, np.where(np.diff(non_zero_index) != 1)[0]+1)
    non_cov_pos_block = np.split(zero_index, np.where(np.diff(zero_index) != 1)[0]+1)
    print (cov_pos_block,non_cov_pos_block)

    # string concatenate

    defalt += ','.join([str(amino_ele_pos_dict[each[0]+1][0])+'-'+str(amino_ele_pos_dict[each[-1]+1][-1])
                      for each in non_cov_pos_block])
    covered += ','.join([str(amino_ele_pos_dict[each[0]+1][0])+'-'+str(amino_ele_pos_dict[each[-1]+1][-1])
                      for each in cov_pos_block])

    # ptm color string concatenate
    ptm_color = ''
    if ptm_idx_dict:
        for ptm in ptm_idx_dict:

            ptm_color += 'color:'+','.join(['%.3f' % (each/256) for each in regex_color_dict[ptm]])+':'
            ptm_color += ','.join([str(amino_ele_pos_dict[idx+1][0]) + '-'
                                   + str(amino_ele_pos_dict[idx+1][-1])
                                   for idx in ptm_idx_dict[ptm]])
            ptm_color += '\n'

    return defalt+'\n'+covered+'\n'+ptm_color.rstrip('\n')


def dump_rep_color_from_array(name, freq_array, ptm_idx_dict, regex_color_dict, base_path=None):
    """
    add color rep without Pymol API
    :param name:
    :param freq_array:
    :param ptm_idx_dict:
    :param regex_color_dict:
    :param base_path:
    :return:
    """
    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'), os.path.join(os.getcwd(), 'js'))
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
            print (ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

    pdb_str = cmd.get_pdbstr(name)
    ret += '\n'+color_getter(freq_array,ptm_idx_dict,regex_color_dict,pdb_str)

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
        template = open(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'imported.html')).read().\
            replace("###INCLUDE_PDB_FILE_HERE###", pdb_str).\
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)
    else:
        template = open('imported.html').read().\
            replace("###INCLUDE_PDB_FILE_HERE###", pdb_str).\
            replace('###INCLUDE_REPRESENTATION_HERE###',ret)

    if base_path:
        f = open(base_path+name + '.html', 'w')
        print (f'html file to {base_path+name}.html')
    else:
        f = open(name.split('-')[1] + '.html', 'w')
        print ('html file to %s' % name.split('-')[1]+'.html')
    f.write(template)
    f.close()


def dump_rep_server(name, freq_array, ptm_idx_dict, regex_color_dict, base_path, bg_color):
    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'), os.path.join(os.getcwd(), 'js'))
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
            print (ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

    pdb_str = cmd.get_pdbstr(name)
    ret += '\n'+color_getter(freq_array,ptm_idx_dict,regex_color_dict,pdb_str)

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

    ret += "\nbgcolor:"+bg_color
    # bgcolor = cmd.get_setting_tuple('bg_rgb')[1]
    #
    # if len(bgcolor) == 1:
    #     bgcolor = cmd.get_color_tuple(bgcolor[0])
    #
    # ret += "\nbgcolor:%02x%02x%02x" % (int(255 * float(bgcolor[0])), \
    #                                    int(255 * float(bgcolor[1])), int(255 * float(bgcolor[2])))
    # if 'PYMOL_GIT_MOD' in os.environ:
    #     template = open(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'imported.html')).read().\
    #         replace("###INCLUDE_PDB_FILE_HERE###", pdb_str).\
    #         replace('###INCLUDE_REPRESENTATION_HERE###', ret)
    # else:
    #     template = open('imported.html').read().\
    #         replace("###INCLUDE_PDB_FILE_HERE###", pdb_str).\
    #         replace('###INCLUDE_REPRESENTATION_HERE###',ret)
    #
    # if base_path:
    #     f = open(base_path+name + '.html', 'w')
    #     print (f'html file to {base_path+name}.html')
    # else:
    #     f = open(name.split('-')[1] + '.html', 'w')
    #     print ('html file to %s' % name.split('-')[1]+'.html')
    # f.write(template)
    # f.close()
    dict = {'pbdstr': cmd.get_pdbstr(name), 'ret': ret}
    with open(base_path + '.json', 'w') as f:
        json.dump(dict, f)


def pdb2json(pdb_file,base_path):
    import pymol
    import json

    pdb_name = os.path.split(pdb_file)[1]
    print(pdb_name)
    pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    pymol.finish_launching()
    pymol.cmd.load(pdb_file, pdb_name)
    # pymol.cmd.disable("all")
    # pymol.cmd.enable()
    # pymol.cmd.hide('all')
    # pymol.cmd.show('cartoon')

    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'], 'pymol2glmol', 'js'), os.path.join(os.getcwd(), 'js'))
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
        if (obj[1] == 0 and obj[4] == 1 and obj[0] == pdb_name):
            ret += parseObjMol(obj)
            # print (ret)
        if (obj[1] == 0 and obj[4] == 4):  # currently all dist objects are exported
            ret += parseDistObj(obj)

    pdb_str = cmd.get_pdbstr(pdb_name)

    view_str = ''
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
    view_str += "\nview:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % \
           (cx, cy, cz, cameraZ, slabNear, slabFar, fogStart, fov)
    for i in range(9):
        view_str += ",%.3f" % view[i]
    dict = {'pdb_str':pdb_str,'ret':ret,'view':view_str}

    with open(base_path+'/'+pdb_name.split('-')[1]+'.json', 'w') as f:
        json.dump(dict,f)


if __name__ == '__main__':
    import numpy as np
    import re
    import time


    freq_array = np.array([0,0,0,1,1,1,1,0,0,0,0,2,2,2,1,1,1,0,0,0])
    # color_str = color_getter(freq_array)
    # print (color_str)
    pdb_file = 'D:/data/alphafold_pdb/UP000000589_10090_MOUSE/AF-Q8TER0-F1-model_v1.pdb'
    time_start = time.time()
    pdb2json(pdb_file,base_path=os.getcwd())