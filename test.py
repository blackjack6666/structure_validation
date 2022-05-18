import pickle
import freesasa
from morfeus import SASA, read_xyz
from pdb_operation import atom_density_center

# chymo_cleav_density_dict = pickle.load(open('D:/data/alphafold_pdb/688_prot_chymotry_cleave_density_dict.pkl', 'rb'))
# print(chymo_cleav_density_dict['P04406'])

# structure = freesasa.Structure("C:/Users/gao lab computer/Downloads/AF-Q96AE4-F1-model_v2.pdb")
# result = freesasa.calc(structure)
# print (result.residueAreas()['A']['106'].total )
#
# elements, coordinates = read_xyz("C:/Users/gao lab computer/Downloads/dodecahedron.xyz")
# print (elements,coordinates)


pick_result = pickle.load(open('D:/data/alphafold_pdb/1207control_protein_KR_sasa_dict.pkl', 'rb'))

for each in pick_result['Q96AE4']:
    print(each, pick_result['Q96AE4'][each])
    # print(each, pick_result['Q96AE4'][each], sasa.atom_areas[each])

# print (atom_density_center('C:/Users/gao lab computer/Downloads/AF-Q96AE4-F1-model_v2.pdb',radius=5.0))
