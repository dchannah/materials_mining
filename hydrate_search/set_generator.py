#!/usr/bin/python

import yaml
import urllib
from pymongo import MongoClient
from pymatgen import MPRester


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


def generate_query(pick1_list, must_have_list, ne=3):
    """
    MongoDB query for any element from pick1 list and all from must have.
    :param pick1_list: A list of elements from which one must be present.
    :param must_have_list: All elements in this list must be present
    :param ne: Number of elements in the compound
    :return: A list of documents meeting those criteria
    """
    e_hull_cutoff = 0.050
    ehq = {"e_above_hull": {"$lt": e_hull_cutoff}}
    s_p1 = str(pick1_list)
    s_mh = str(must_have_list)
    ne = str(ne)
    # cq = "{\"elements\": {\"$in\": " + s_p1 + ", \"$all\": " + s_mh + \
    #      "}, \"nelements\": " + ne + "}"
    cq = "{\"elements\": {\"$in\": " + s_p1 + ", \"$all\": " + s_mh + \
         "}}"
    tq = "{ \"$and\": [" + str(ehq) + ", " + cq + "]}"
    tq = eval(tq)
    print("ABOUT TO QUERY " + str(tq))
    with MPRester(api_key="cVkBCOZrv4TehHfw") as m:
        data = m.query(tq, ["material_id", "pretty_formula"])
        # data = m.query(tq, ["snl_final"])
    return data


# yaml_txt = "/Users/dchannah/google_drive/yaml/empty_local.yaml"
# mc = initialize_db(yaml_txt)
# leach_metals = ["Ag"]
# anion = ["O"]
# tms = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Mo"]
# q = generate_query(tms, leach_metals + anion, 4)

# for lm in leach_metals:
#     lm_l = [lm]
#     print "Generating a query for", lm_l, "as a leachable ion."
#     # q = generate_query(tms, lm_l + anion, 3)
#     q = generate_query(tms, anion, 2)
#     for cmpd in q:
#         mc.insert_one(cmpd)
#     print "Inserted compounds for", lm, "successfully."
