import Bio.PDB
import numpy
from Bio.PDB import PDBIO

def residue_dist_to_ligand(protein_residue, ligand_residue) :
    #Returns distance from the protein C-alpha to the closest ligand atom
    dist = []
    for atom in ligand_residue :
        if "CA" in protein_residue:
            vector  = protein_residue["CA"].coord - atom.coord
            dist.append(numpy.sqrt(numpy.sum(vector * vector)))
            return min(dist)

def get_ligand_by_name(residue_name, model):
    #Extract ligands from all chains in a model by its name
    global ligands
    ligands = {}
    chains = model.child_dict
    for c in chains:
        ligands[c] = []
        for protein_res in chains[c].child_list:
            if protein_res.resname == residue_name:
                ligands[c].append(protein_res)
    return ligands

def get_ligand_by_chain(chain_name, model):
    #Extract all ligand residues from given chain name
    global ligands
    ligands = {}    
    ligands[chain_name] = []
    chains = model.child_dict
    for protein_res in chains[chain_name].child_list:
        ligands[chain_name].append(protein_res)
    return ligands

def active_site(ligands, distance, model):
    # Prints out residues located at a given distance from ligand
    chains = model.child_dict
    for group in ligands.values():
        for ligand_res in group:
            print ligand_res.resname, ligand_res.id[1]
            for c in chains:
                for protein_res in chains[c].child_list:
                    if protein_res not in group:
                        dist = residue_dist_to_ligand(protein_res, ligand_res)
                        if dist and dist < distance :
                            print protein_res.resname, protein_res.id[1], dist


def save_ligand(structure, filename):
    # Saves ligand to a filename.pdb
    Select = Bio.PDB.Select
    class LigandSelect(Select):
        def accept_residue(self, residue):
            for group in ligands.values():
                if residue in group:
                    return 1
                else:
                    return 0
    io=PDBIO()
    io.set_structure(structure)
    io.save(filename+'.pdb', LigandSelect())

"""
Example
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
"""
