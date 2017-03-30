import yaml, urllib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import *
from pymongo import MongoClient
from database_getter import call_derivation, stored, derived
from scipy.stats import linregress

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
    # Initialize the paper2 ionic mobility database.
    db_conf = "/Users/dchannah/yaml/paper2.yaml"
    coll = initialize_db(db_conf)

    # Get all the compounds from the database containing magnesium.
    db_query = build_db_query("cation_type", "Ca")
    # db_query = None
    cathode_documents = coll.find(db_query)

    # Define a list of properties we wish to examine.
    whole_path_props = ['min_anion_distance', 'volume_per_anion', 
                        'min_cation_distance', 'site_cn_diff', 
                        'average_bond_length', 'av_delta_cn']

    fv_x = build_feature_vector(whole_path_props, cathode_documents)
    
    # Get the value we wish to regress against the properties above.
    cathode_documents = coll.find(db_query)
    y_prop = ['activation_energy']
    fv_y = build_feature_vector(y_prop, cathode_documents)

    # Create subplots showing the x vs. y data for each property.
    fig = plt.figure(figsize=(146, 8))
    num_props = len(whole_path_props)
    cols = num_props/2
    rows = int(math.ceil(num_props/cols))
    gs = gridspec.GridSpec(rows, cols)
    
    # Create a figure to show how well the fits are doing.
    fig_err = plt.figure()
    bar_labels = []
    errors = []
    for i in range(num_props):
        # Need to get data in xy format.
        x_data, y_data = [], []
        for path in fv_x:
            x_data.append(fv_x[path][whole_path_props[i]])
            y_data.append(fv_y[path][y_prop[0]])

        ax = fig.add_subplot(gs[i]) 
        ax.scatter(x_data, y_data)

        # Linear fit to look any obvious correlations.
        slope, intercept, r_value, p_value, std_err = linregress(x_data, y_data)
        fit_line = [slope*x + intercept for x in x_data]
        ax.plot(x_data, fit_line)  # Plotting the fit line 

        # Catalog our errors to see which property is best.
        error_metric = r_value**2
        # error_metric = std_err
        bar_labels.append(whole_path_props[i])
        errors.append(error_metric)

        # We need some labels.
        ax.set_xlabel(whole_path_props[i])
        ax.set_ylabel(y_prop[0])

    
    # In a separate figure, we plot the errors.
    ax_error = fig_err.add_subplot(111)
    xtick_positions = [pos + 0.5 for pos in range(num_props)]
    ax_error.set_xticks(xtick_positions)
    ax_error.set_xticklabels(bar_labels, rotation=270)
    ax_error.bar(range(num_props), errors)

    # Show the plot
    plt.tight_layout()
    plt.show()

