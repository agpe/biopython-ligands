biopython-ligands
=================
With this module you can extract ligand from your pdb file and see which residues from receptor are in close relation with it. Also, extracted ligand can be saved to a new pdb file. This module is using BioPDB and its Structure object architecture (the so-called SMCRA (Structure/Model/Chain/Residue/Atom))  :
A structure consists of models
A model consists of chains
A chain consists of residues
A residue consists of atoms  

Functions:
residue_dist_to_ligand(protein_residue, ligand_residue)
Protein_residue and ligand_residue are Residue Objects from BioPDB module.

get_ligand_by_name(residue_name, model)
Residue_name is a string containing name of wanted ligand, model is a Model Object. This function returns a dictionary which is necessary in active_site and save_ligand functions 

get_ligand_by_chain(chain_name, model)
Chain_name is a string containing name of a chain, model is a Model Object. This function returns a dictionary which is necessary in active_site and save_ligand functions.

active_site(ligands, distance, model)
Ligands is a dictionary where keys are chains and values are ligand residues located in that chain, distance is a maximum distance you want to obtain receptor residues from, model is a Model Object.

save_ligand(structure, filename)
Structure is a Structure Object, filename is a name of a pdb file that ligand will be saved to
Example:

pdb_code = "1E8W"
pdb_filename = "1E8W.pdb"
structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
model = structure[0]
ligand = get_ligand_by_name('QUE', model)
active_site(ligand, 5.0, model)

structure2 = Bio.PDB.PDBParser().get_structure('4GRV', '4GRV.pdb')
model2 = structure2[0]
ligand2 = get_ligand_by_chain('B', model2)
active_site(ligand2, 8.0, model2)
save_ligand(structure2, 'ligand')
