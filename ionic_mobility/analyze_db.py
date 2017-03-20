import yaml, urllib
from pymongo import MongoClient
from database_getter import call_derivation, stored, derived


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


def build_feature_vector(prop_list, docs):
    """build_feature_vector
    Constructs a feature vector using a list of properties for a list of docs.
    :param prop_list: List of properties (must be defined in DB or derived)
    :param docs: List of MongoDB documents.
    :return: Dictionary of {"path_id": feature vector}
    """
    # Note that path_ids are UNIQUELY defined in our MongoDB.
    # We don't check this constraint in the script, so be careful.
    all_feature_vectors = {}
    for doc in docs:
        feature_vector = {}
        for prop in prop_list:
            if prop in stored:
                feature_vector[prop] = doc["NEB_analysis"][prop]
            elif prop in derived:
                feature_vector[prop] = call_derivation(prop, doc)
            else:
                feature_vector[prop] = doc[prop]
        all_feature_vectors[doc["path_id"]] = feature_vector
    return all_feature_vectors


def get_traj_props(traj_prop, docs):
    """get_traj_props
    Returns trajectory-resolved properties for a list of documents.
    :param traj_prop: A property defined for each image in the NEB trajectory.
    :param docs: A list of MongoDB documents.
    :return: Dictionary {"path_id": [traj_prop_img_0, ..., traj_prop_img_N]}
    """
    # Need to add-in a check on whether a trajectory-resolved property is
    # actually being requested.
    traj_prop_dict = {}
    for doc in docs:
        if traj_prop in stored:
            traj_prop_dict[doc["path_id"]] = doc["NEB_analysis"][traj_prop]
        elif traj_prop in derived:
            traj_prop_dict[doc["path_id"]] = call_derivation(traj_prop, doc)
        else:
            traj_prop_dict[doc["path_id"]] = doc[traj_prop]
    return traj_prop_dict


def get_doc_props(doc_prop, docs):
    """get_doc_props
    Returns value of a document-resolved value for a list of documents.
    :param doc_prop: A property defined uniquely for a whole NEB trajectory. 
    :param docs: A list of MongoDB documents.
    :return: Dictionary {"path_id": doc_prop}
    """
    # Need to add in a check on whether a doc-resolved prop is requested.
    doc_prop_dict = {}
    for doc in docs:
        if doc_prop in stored:
            doc_prop_dict[doc["path_id"]] = doc["NEB_analysis"][doc_prop]
        elif doc_prop in derived:
            doc_prop_dict[doc["path_id"]] = call_derivation(prop, doc)
        else:
            doc_prop_dict[doc["path_id"]] = doc[doc_prop]
    return doc_prop_dict


def build_db_query(criterion=None, value=None, ucq=None, custom=False):
    """build_db_query
    Generates a MongoDB query for paticular NEB runs.
    :param criterion: Search criterion.
    :param value: Criterion value to search for.
    :param ucq: User-defined custom query.
    :param custom: Boolean - is user defined query being passed?
    """
    if criterion == None and value == None and ucq == None:
        print("You didn't query anything!")
        return
    elif not custom:
        q_s = "{\"" + criterion + "\": \"" + value + "\"}"
        return eval(q_s)  # Eval is a little sketchy, but it'll do for now.
    else:
        return eval(ucq)

if __name__ == "__main__":
    # A simple set of commands to test the script.
    # Let's compare the minimum cation distance for the post-spinels.
    path_id_cm = "mvc-6921-NEB-001"
    path_id_cf = "mvc-12110-NEB-001"
    cm_query = build_db_query("path_id", path_id_cm)
    cf_query = build_db_query("path_id", path_id_cf)
    yaml_file = "/Users/dchannah/yaml/paper2.yaml"
    coll = initialize_db(yaml_file)
    cm_doc = coll.find(cm_query)
    cf_doc = coll.find(cf_query)
    print "Marokite:", get_traj_props("min_cat_dist", cm_doc)
    print "CaFe2O4:", get_traj_props("min_cat_dist", cf_doc)
