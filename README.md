# Alpha-VAL
## Scripts for alphafold2 structure validation from native digestion
### Some functions for computation

### 1. solvent accessibility surface area 
Calculated by PyMOL API pymol.cmd.get_area
Implemented in pdb_operation.py, "sasa_pdb", line 456-484, executed in line 825-837
```python
import multiprocessing

start = time.time()
with multiprocessing.Pool(multiprocessing.cpu_count() - 2) as pool:
    result = pool.map(sasa_pdb, input_list_tuples, chunksize=50)
    pool.close()
    pool.join()
file_density_dict = {k: v for d in result for k, v in d.items()}

pickle.dump(file_density_dict, open('D:/data/alphafold_pdb/1207control_protein_KR_sasa_dict.pkl', 'wb'))
print(time.time() - start)
```

### 2. atom density (number of proximal atoms near cleavage site)
find number of atoms in a sphere space for each in-silico digested site
Code: pdb_operation.py, functions "residue_density_cal2", "inSphere2", executed in line 805

### 3. cleavage to geometric center distance
calculate the euclidean distance between cleavage and geometric center of protein
Code: pdb_operation.py, "find_centroid", "residue_distance", calculation in batch: line 595
```python
def residue_distance(residue_atom_xyz):
    """
    calculate the distance for each residue to center of 3d structure
    :param residue_atom_xyz:
    :return:
    """
    # zero_point = np.array([0,0,0])
    zero_point = find_centroid(residue_atom_xyz) # function to find genometic center
    residue_distance_dict = {}
    for each_pos in residue_atom_xyz:

        total_dist = sum([np.linalg.norm(np.array(each_atom)-zero_point)
                          for each_atom in residue_atom_xyz[each_pos]]) # total euclidean distance of one residue

        average_dist = total_dist/len(residue_atom_xyz[each_pos]) # average by number of atoms
        residue_distance_dict[each_pos] = average_dist
    return residue_distance_dict
```
### other analysis scripts to output graphs are in project "deep_proteome", "structure_validation.py" and "native_digestion_data_analysis.py"
### Contributing
### Authors and acknowledgment
### License