from Bio.PDB import PDBList, PDBParser
import pandas as pd
import numpy as np
from Bio.PDB.vectors import calc_dihedral


def calculate_torsion_angle(atom1, atom2, atom3, atom4):
    vector1 = atom1.get_vector()
    vector2 = atom2.get_vector()
    vector3 = atom3.get_vector()
    vector4 = atom4.get_vector()
    angle = calc_dihedral(vector1, vector2, vector3, vector4)
    angle_in_degrees = angle * (180.0 / np.pi)
    return angle_in_degrees


def extract_torsion_angles(structure):
    torsion_angles_data = []

    for model in structure:
        for chain in model:
            # filtrowanie heteroatomów
            filtered_chain = [residue for residue in chain if residue.id[0] == ' ']
            for i, residue in enumerate(filtered_chain):
                if residue.get_resname() in ['A', 'U', 'G', 'C']:
                    angles = {}
                    try:
                        if i > 0:  # kąt alfa, porzednia reszta
                            prev_residue = filtered_chain[i - 1]
                            angles['alpha'] = calculate_torsion_angle(prev_residue["O3'"], residue["P"], residue["O5'"],
                                                                      residue["C5'"])

                        # pozostałe kąty
                        angles['beta'] = calculate_torsion_angle(residue["P"], residue["O5'"], residue["C5'"],
                                                                 residue["C4'"])
                        angles['gamma'] = calculate_torsion_angle(residue["O5'"], residue["C5'"], residue["C4'"],
                                                                  residue["C3'"])
                        angles['delta'] = calculate_torsion_angle(residue["C5'"], residue["C4'"], residue["C3'"],
                                                                  residue["O3'"])
                        if i < len(filtered_chain) - 1:
                            next_residue = filtered_chain[i + 1]
                            angles['epsilon'] = calculate_torsion_angle(residue["C4'"], residue["C3'"], residue["O3'"],
                                                                        next_residue["P"])
                            angles['zeta'] = calculate_torsion_angle(residue["C3'"], residue["O3'"], next_residue["P"],
                                                                     next_residue["O5'"])

                        # kąty dla cukru
                        angles['v0'] = calculate_torsion_angle(residue["C4'"], residue["O4'"], residue["C1'"],
                                                                residue["C2'"])
                        angles['v1'] = calculate_torsion_angle(residue["O4'"], residue["C1'"], residue["C2'"],
                                                                residue["C3'"])
                        angles['v2'] = calculate_torsion_angle(residue["C1'"], residue["C2'"], residue["C3'"],
                                                                residue["C4'"])
                        angles['v3'] = calculate_torsion_angle(residue["C2'"], residue["C3'"], residue["C4'"],
                                                                residue["O4'"])
                        angles['v4'] = calculate_torsion_angle(residue["C3'"], residue["C4'"], residue["O4'"],
                                                                residue["C1'"])

                        # kąt chi dla zasady azotowej
                        base_atoms = {
                            'A': ('N9', 'C4'),
                            'G': ('N9', 'C4'),
                            'C': ('N1', 'C2'),
                            'U': ('N1', 'C2')
                        }
                        base_atom, prev_atom = base_atoms[residue.get_resname()]
                        angles['chi'] = calculate_torsion_angle(residue["O4'"], residue["C1'"], residue[base_atom],
                                                                residue[prev_atom])

                        # kąty do listy
                        torsion_angles_data.append(angles)
                    except KeyError as e:
                        print(f"Residue {residue.id[1]} is missing an atom for torsion calculation: {e}")
                        continue

    return pd.DataFrame(torsion_angles_data)


def load_structure(pdb_id):
    pdb_list = PDBList()
    pdb_path = pdb_list.retrieve_pdb_file(pdb_id, pdir='.', file_format="pdb")
    parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, pdb_path)


# pdbidcząsteczki
pdb_id = "1EHZ"

# wczytanie
structure = load_structure(pdb_id)

# ekstrakcja kątów torsyjnych
torsion_angles_df = extract_torsion_angles(structure)

# zapis
csv_path = f'{pdb_id}_torsion_angles.csv'

# pzekształcenie danych do postaci macierzy
torsion_angles_matrix = torsion_angles_df.values.reshape(-1, len(torsion_angles_df.columns))
columns_names = list(torsion_angles_df.columns)

# zamień NaN na zero
# torsion_angles_matrix[np.isnan(torsion_angles_matrix)] = 0

# zapis
np.savetxt(csv_path, torsion_angles_matrix, delimiter=",", fmt="%.3f", header=",".join(columns_names), comments="")
