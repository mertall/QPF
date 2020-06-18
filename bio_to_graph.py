from Bio.PDB import PDBParser
import networkx as nx
import numpy as np
import openbabel
import sys
import matplotlib.pyplot as plt

# Recommended pipelines:
#	2D sdf --create_pdb--> 2D pdb --create_structure--> Biopython structure --create_graph--> 2D networkx graph
#	3D sdf/pdb --create_pdb--> 2D pdb without hydrogens --create_structure--> Biopython structure --create_graph--> 2D networkx graph, no H's
# 	3D pdb --antechamber, not included here--> mpdb with charges --create_structure--> Biopython structure --create_graph--> 3D networkx graph


# Convert 3D pdb to 2D pdb, while removing hydrogens because openbabel doesn't support hydrogens in 2D drawings
def twoD_eliminateH(filename):
	pybel.readfile("pdb", filename)
	newfilename = "twoD_eliminateH_"+filename
	output = pybel.Outputfile("pdb", newfilename)
	for mol in pybel.readfile("pdb", filename):
		mol.removeh()
		mol.draw(show = False, update = True)
		output.write(mol)
	return newfilename

# FUNCTION: create_pdb
# ARGUMENTS: filename (string with extension included), format (string), quiet (bool)
# 	filename
#	a_format. Accepted formats: any accepted by openBabel's OBConversion. "smi" input should have no hydrogens. Only smi, pdb, and sdf tested.
# 	create_twoD_eliminateH. Set to True only for 3d-->2d conversion. Result will exclude hydrogens
# OUTPUT: a PDB file is written

def create_pdb(filename, a_format, create_twoD_eliminateH = False):
	
	fileid = filename.rsplit(".", 1)[0]
	print("checking filetype and performing conversions...")

	# Convert filetypes to to pdb, write pdb file and use that file for rest of function
	if a_format != "pdb":
		try:
			obConversion = openbabel.OBConversion()
			obConversion.SetInAndOutFormats(a_format, "pdb")
			mol = openbabel.OBMol()
			obConversion.ReadFile(mol, filename)   # Open Babel will uncompress automatically if given .gz
		except:
			print("Open babel could not read file for conversion. Check filename, format, and file content")

		if a_format == "smi":
			mol.AddHydrogens()

		obConversion.WriteFile(mol, fileid+".pdb")
		filename = fileid+".pdb"

	# Reduce to 2-dimensional. Will write another pdb file and use that file for rest of function
	if create_twoD_eliminateH:
		print("creating twoD_eliminateH pdb")
		filename = twoD_eliminateH(filename)

	print("pdb file created: " + filename)
	return filename

# FUNCTION: create_structure
# ARGUMENTS: pdb filename. If you want to include charges, you have to add them yourself into the occupancy
# column (the column to the right of the z coordinate). 
# OUTPUT: Biopython structure 
def create_structure(filename, quiet = True):
	print("creating biopython molecule structure...")
	fileid = filename.rsplit(".", 1)[0]
	p = PDBParser(QUIET=quiet) 
	structure = p.get_structure(fileid, filename)
	return structure

def bond_length(atom1, atom2):
	# bond lengths from LibreText (added 0.1 for variability) https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Chemical_Bonds/Bond_Lengths_and_Energies
	# units: angstroms. 
	# TODO: make bond lengths include sulfur. Make more efficient (change data structure?)
	# bond lengths for same atoms (angstroms)
    if atom1.element == "C" and atom2.element=="C":
            return 1.64
    elif atom1.element == "O" and atom2.element=="O":
            return 1.58
    elif atom1.element == "N" and atom2.element=="N":
            return 1.3
    elif atom1.element == "H" and atom2.element=="H":
            return 0.84

    #bond lengths for all possible atom pairs with C, O, N, and H. Users will input for other elements. Most amino acids will just have this.
    elif atom1.element == "C" and atom2.element=="O":
            return 1.53
    elif atom1.element == "O" and atom2.element=="C":
            return 1.53
    elif atom1.element == "C" and atom2.element=="H":
            return 1.18
    elif atom1.element == "H" and atom2.element=="C":
            return 1.18
    elif atom1.element == "C" and atom2.element=="N":
            return 1.57
    elif atom1.element == "N" and atom2.element=="C":
            return 1.57
    elif atom1.element == "H" and atom2.element=="N":
            return 1.11
    elif atom1.element == "N" and atom2.element=="H":
            return 1.11
    elif atom1.element == "O" and atom2.element=="N":
            return 1.3
    elif atom1.element == "N" and atom2.element=="O":
            return 1.3
    elif atom1.element == "O" and atom2.element=="H":
            return 1.06
    elif atom1.element == "H" and atom2.element=="O":
            return 1.06
    else:
            return float(input("Insert an estimate bond length for " + str(atom1.element) + " to " + str(atom2.element) + " in Angstroms: "))
        
# FUNCTION: create_graph
# ARGUMENTS: structure (biopython class), has_charge is True unless you are starting from a pdb (not mpdb) or sdf
# OUTPUT: networkx graph. 
# NOTES: Node items are atoms with attributes charge and pos. Fully connected graph. Edges have attribute len. 
def create_graph(structure, has_charge=True):
	# Graph atoms
	G = nx.Graph()
	locations = []
	names = []
	# Initialize nodes with atom name, charge, and position
	print("creating graph nodes...")
	atoms = structure.get_atoms()
	for atom in atoms:
			location = atom.get_coord() #converted to string for graphml export
			locations.append(location)
			names.append(atom.get_name())

			if has_charge:
				bcc_charge = atom.get_occupancy()
				G.add_node(atom.get_name(), charge=bcc_charge, x=location[0], y=location[1], z=location[2]) # could also add pos=location, doesn't work for export?
			else:
				G.add_node(atom.get_name(), charge=0, x=location[0], y=location[1], z=location[2])

	# Create edges based on distance. Calculate distances for all atom pairs.
	# TODO: When scaling up, change atom2 loop to only be atoms within a certain radius of atom1. For now this will do.
	
	print("creating graph edges...")


		
		

	for atom1 in structure.get_atoms():
		for atom2 in structure.get_atoms():
			distance = atom1 - atom2
			if (0 < distance <= bond_length(atom1, atom2)) and not ((atom1.element and atom2.element) == "H"):
				G.add_edge(atom1.get_name(), atom2.get_name(), len=distance) # subtraction operator calculates distance. {'weight': atom1-atom2}

	return G,locations,names

def visualize_graph(graph):
	nx.draw(graph, with_labels=True)
	plt.show()

def main():
	'''
	# Create graph from sdf, export as graphml
	filename = "gly.sdf"
	filename = create_pdb(filename, "sdf")
	structure = create_structure(filename)
	graph = create_graph(structure, has_charge = False)
	nx.write_graphml(graph, "gly.graphml") # change pathname for your environment
	copy = nx.read_graphml("gly.graphml")
	print("original adjacency representation")
	print(graph.nodes(data=True))
	print(graph.edges(data=True))
	print("\ncopied adjacency representation")
	print(copy.nodes(data=True))
	print(copy.edges(data=True))
	visualize_graph(graph)
	'''
if __name__ == '__main__':
	main()