#!/usr/bin/python

from new_db_field import initialize_db, id_mobile_ion
from pymatgen import Structure, Element, MPRester, Specie
from pymatgen.analysis.structure_matcher import StructureMatcher

# Global variables

dan_yaml = "/Users/dchannah/yaml/empty_local.yaml"
base_db_yaml = "/Users/dchannah/yaml/pf_db.yaml"
mobile_ion = ["Li", "Na", "Mg", "Ca", "Zn", "K", "Ag", "Cu"]
api = "cVkBCOZrv4TehHfw"

# Initialize the collections
dan_coll = initialize_db(dan_yaml)
base_db_coll = initialize_db(base_db_yaml)

# Get a dict of all the gamma structures from Miao's collection.
base_db_all_docs = base_db_coll.find()
gamma_structures = {}
for doc in base_db_all_docs:
    gamma_structures[doc["mp-id"]] = Structure.from_dict(doc["gamma_structure"])

# Get all the structures in Dan's collection
dan_all_docs = dan_coll.find()
dan_empty_structures = {}
for doc in dan_all_docs:
    with MPRester(api) as m:
        s_doc = m.get_structure_by_material_id(doc["material_id"])
    # formula_str = str(s_doc.composition.reduced_composition).replace(" ", "")
    # working_ion = id_mobile_ion(formula_str)
    # s_doc.remove_species([working_ion])
    dan_empty_structures[doc["material_id"]] = s_doc

# Check for matches
print "Starting structure matching"
sm = StructureMatcher(primitive_cell=False, scale=True, attempt_supercell=True)
for mpid, struct in dan_empty_structures.iteritems():
    for base_db_mpid, base_db_struct in gamma_structures.iteritems():
        if sm.fit(struct, base_db_struct):
            print mpid, "matches", base_db_mpid
