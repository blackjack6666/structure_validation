import pickle

chymo_cleav_density_dict = pickle.load(open('D:/data/alphafold_pdb/688_prot_chymotry_cleave_density_dict.pkl', 'rb'))
print(chymo_cleav_density_dict['P04406'])
