biopython-ligands
=================
With this module you can extract ligand from your pdb file and see which residues from receptor are in close relation with it. Also, extracted ligand can be saved to a new pdb file.   


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
