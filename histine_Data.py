import Bio.PDB
import numpy as np
import sys
from bio_to_graph import create_pdb,create_structure, create_graph # create_pdb(), create_structure(), and create_graph()


def grab_db():
    filename = "his.sdf"
    filename = create_pdb(filename, "sdf")
    structure = create_structure(filename)
    graph,locations,names = create_graph(structure, has_charge =True)
    print(graph.edges(data=True))
    unformatted_DB = locations
    print(str(locations))
    print(type(unformatted_DB))
    return unformatted_DB,names
print(grab_db())
