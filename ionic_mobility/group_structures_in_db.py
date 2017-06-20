#!/usr/bin/python

import os, yaml, sys, urllib, ast
from pymongo import MongoClient
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.structure_matcher import FrameworkComparator


def initialize_db(yaml_file):
    """
    Initializes a MongoClient object with a user-specified database.
    :param yaml_file: YAML file containing database configuration.
    :return: A Mongo object with the initialized database.
    """
    with open(yaml_file, 'r') as f:
        db_config = yaml.load(f)
    client = MongoClient(
        "mongodb://" + db_config['user'] + ":" + urllib.quote_plus(
            db_config['passwd']) + "@" + db_config['host'] + ":" + str(
            db_config['port']) + "/" + db_config['db'])
    db = client[db_config['db']]
    collection = db[db_config['collection']]
    return collection


def query_db(criterion, value, coll):
    """query_db
    Queries a MongoDB and returns a list of matching documents.
    :param criterion: Key for the search dict.
    :param value: Value for the search dict.
    :param coll: A Pymongo collection.
    :return: A List of matching PyMongo Documents.
    """
    query = "{\"" + criterion + "\": \"" + value + "\"}"
    query = ast.literal_eval(query)
    return coll.find(query)


def get_empty_structures(criterion, s_values, coll):
    """get_empty_structures
    Gets all host structures associated with a list of search values.
    :param criterion: The search criterion to use to get structures.
    :param s_values: List of search values (e.g. path IDs) to get data for.
    :param coll: A PyMongo collection.
    :return: A dictionary of {s_value: empty_structure} 
    """
    s_dict = {}
    for s in s_values:
        results = query_db(criterion, s, coll)
        structure = Structure.from_dict(results[0]["host_structure"])
        s_dict[s] = structure
    return s_dict


def list_from_dict(d):
    """list_from_dict
    Converts a dict to a list - useful for data processing.
    :param d: A dictionary.
    :return: Lists of keys and values with same indexing.
    """
    l_k, l_v = [], []
    for key in d:
        l_k.append(key)
        l_v.append(d[key])
    return l_k, l_v


def group_structures(s_l):
    """group_structures
    Applies a structure grouping algorithm to a list of structures.
    :param s_l: List of Structure objects.
    :return: List of lists of grouped structures.
    """
    sm = StructureMatcher(scale=True, attempt_supercell=True, 
                          comparator=FrameworkComparator())
    return sm.group_structures(s_l)


def read_list(f_path):
    """read_list
    Converts lines of a file into a list.
    :param f_path: Path to file.
    :return: List of the lines in the file.
    """
    l = []
    with open(f_path, 'r') as f:
        for line in f:
            l.append(line.rstrip())
    return l


def get_ids_in_bin(grouped_l, s_dict):
    """get_ids_in_bin
    Gets a list of search values (path id, for example) in each struct bin.
    :param grouped_l: A list of lists of grouped structures.
    :param s_dict: A structure dict mapping a struct to a search val.
    :return: A dictionary of bins and the search values in that bin.
    """
    structure_bins = {}
    num_groups = len(grouped_l)
    for i in range(num_groups):
        structure_bins[i] = []
    for group_number, sub_group in enumerate(grouped_l):
        for s in sub_group:
            for search_string, struct in s_dict.items():
                if s == struct:
                    structure_bins[group_number].append(search_string)
    return structure_bins


yaml_file = "/Users/dchannah/yaml/paper2.yaml"
db = initialize_db(yaml_file)
path_id_list_file = "./search.txt"
path_ids = read_list(path_id_list_file) 
structures = get_empty_structures("path_id", path_ids, db)
searches, structures_list = list_from_dict(structures)
grouped_structures = group_structures(structures_list)
structs_in_groups = get_ids_in_bin(grouped_structures, structures)
for k in structs_in_groups:
    print "Number of structures in group:", len(structs_in_groups[k]), ", IDs in group:", structs_in_groups[k]
