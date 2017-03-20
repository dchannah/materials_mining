#!/usr/bin/python
from diffusion_updater_v3 import UpdaterNEB
import os
from pymongo import MongoClient
from pymatgen import MPRester
import yaml
import sys
import urllib
import time
import numpy as np
import argparse
__author__ = 'dchannah'


def initialize_database(mp_id):
    with open(db_config_path, 'r') as f:
        db_config = yaml.load(f)
    client = MongoClient("mongodb://" + db_config['user'] + ":" + urllib.quote_plus(db_config['passwd']) + "@" +
                         db_config['host'] + ":" + str(db_config['port']) + "/" + db_config['db'])
    db = client[db_config['db']]
    collection = db[db_config['collection']]
    user_info = {"name": db_config['name'], "email": db_config['email']}
    if mp_id[0:6] == "NoMPID":
        pretty_formula = mp_id[6:]
    else:
        with MPRester(api_key=db_config['API']) as m:
            pretty_formula = m.get_data(mp_id, prop='pretty_formula')[0]['pretty_formula']
    return collection, pretty_formula, user_info


def upload_to_database(path_type, mp_id, working_ion, pair):
    mp_collection, pretty_formula, user_info = initialize_database(mp_id)

    if path_type == "NEB":
        update = UpdaterNEB(mp_id, working_ion, pair)
        db_dict = update.populate_fields()
        mp_collection.update_one({"path_id": db_dict['path_id']},
                {"$set": {"path_id": db_dict['path_id'], 
                          "dir": os.getcwd(),
                          "info": {"generated_by": user_info['name'], "email": user_info['email'],
                                   "time": time.strftime("%c")},
                          "mp_id": mp_id,
                          "pretty_formula": pretty_formula,
                          "path_type": db_dict['path_type'],
                          "cation_type": working_ion,
                          "host_structure": db_dict['host_structure'],
                          "cation_concentration": db_dict['cation_concentration'],
                          "is_vol_relaxed": db_dict['vol_relaxed'],
                          "cation_diffusion_start": db_dict['cation_diffusion_start'],
                          "cation_diffusion_end": db_dict['cation_diffusion_end'],
                          "inputs": db_dict['input_set'],
                          "path_sites": db_dict['path_sites'],
                          "neb_images": db_dict['neb_images'],
                          "NEB_analysis": {"path_energy": db_dict['path_energy'],
                                           "path_length": db_dict['path_length'],
                                           "coordination": db_dict['coordination'],
                                           "effective_coordination": db_dict['effective_coordination'],
                                           "average_bond_length": db_dict['average_bond_length'],
                                           "distortion_index": db_dict['distortion_index'],
                                           "activation_energy": db_dict['activation_energy']},
                          "chgcar_status": "success",
                          "valid": "True"}}, True)

    elif path_type == "approxNEB":
        print "Sorry, approxNEB isn't implemented yet."
        # update = updater_approx_neb(mp_id, working_ion)
        # db_dict = update.populate_fields()

    elif path_type == "voronoi":
        print "Sorry, voronoi isn't implemented yet."
        # update = updater_voronoi(mp_id, working_ion)
        # db_dict = update.populate_fields()

    else:
        print "Exiting; path type must be NEB, approxNEB, or voronoi."
        sys.exit()

    return


parser = argparse.ArgumentParser()
parser.add_argument("mp_id", help="Materials Project ID or NoMPID<CompoundName>")
parser.add_argument("calc_type", help="Type of calculation (NEB, approxNEB, or Voronoi)")
parser.add_argument("cation", help="Intercalating cation")
parser.add_argument("pair_index", help="The index of the end-point pair you are calculating.  Internal bookkeeping.")
parser.add_argument("db_config_path", help="Global path to a YAML file containing user info and database config.")

args = parser.parse_args()
id = str(args.mp_id)
calc_type = str(args.calc_type)
cation = str(args.cation)
pair_index = str(args.pair_index)
db_config_path = str(args.db_config_path)

upload_to_database(calc_type, id, cation, pair_index)
