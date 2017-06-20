#!/usr/bin/python

import yaml
import urllib
from oxi_filter import compute_grav_cap
from pymongo import MongoClient
from pymatgen import Composition

f_yaml = "/Users/dchannah/google_drive/yaml/empty_local.yaml"
mobile_ion = ["Li", "Na", "Mg", "Ca", "Zn", "K", "Ag", "Cu"]

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


def update_field(coll, s_c, s_v, f_n, f_v):
    """update_field
    Updates the database given a particular criterion, name, and value.
    :param coll: A PyMongo Collection.
    :param s_c: Search criteria, should return *something*
    :param s_v: Thing to search for.
    :param f_n: Name of the field to create or update.
    :param f_v: Value for the field we are creating/updating.
    :return: None
    """
    coll.update_one({s_c: s_v}, {"$set": {f_n: f_v}}, True)
    return


def get_all_ids(coll):
    """get_all_ids
    Gets all the MPID values from a particular MongoDB collection.
    :param coll: A MongoDB collection
    :return: A list of MPIDs in the collection.
    """
    mpid_list = []
    all_docs = coll.find()
    for d in all_docs:
        mpid_list.append(d["material_id"])
    return mpid_list


def id_mobile_ion(formula):
    """id_mobile_ion
    Given a formula, identifies the mobile ion present, if any.
    :param formula: Formula to fetch working ion from.
    """
    c = Composition(formula)
    wi = [el for el in c.as_dict() if el in mobile_ion]
    return wi[0]


def cap_from_id(coll, mpid, wi=None):
    """cap_from_id
    Given an MPID, fetches the formula and computes capacity from it.
    :param coll: Mongo collection.
    :param mpid: An MPID to get the capacity values for.
    :param wi: String identity of the working ion.
    :return: A gravimetric capacity in mAh/g.
    """
    document = coll.find({"material_id": mpid})[0]
    f = document["pretty_formula"]
    remove = False
    if wi is None:
        wi = id_mobile_ion(f)
        remove = True
    grav_cap = compute_grav_cap(Composition(f), wi, remove)
    return grav_cap

if __name__ == "__main__":
    collection = initialize_db(f_yaml)
    mpids = get_all_ids(collection)
    search_criterion = "material_id"
    f_name = "screening_set_status"
    f_value = "failed"
    
    grav_cap_dict = {}
    for m in mpids:
        grav_cap_dict[m] = cap_from_id(collection, m, "Mg")
        # update_field(collection, search_criterion, m, field_name, field_value)
    
    for mpid in grav_cap_dict:
        if grav_cap_dict[mpid] < 50:
            update_field(collection, search_criterion, mpid, f_name, f_value)
            print "Updated collection for", mpid
            fn2 = "reason_failed"
            fv2 = "Capacity lower than 50 mAh/g"
            update_field(collection, search_criterion, mpid, fn2, fv2)
