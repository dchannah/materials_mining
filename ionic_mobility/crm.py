import numpy as np
import scipy.cluster.vq as vq
from pylab import plot, show
from analyze_db import initialize_db, build_db_query, build_feature_vector

"""
Clustering-ranking-modeling to learn ionic mobility predictors.  We're using
the X-means clustering algorithm here.
"""


def build_point_list(fv_dict):
    """build_point_list
    Builds a list of feature vectors for clustering analysis.
    :param fv_dict: Dictionary of all needed feature vectors. 
    :return: [feature vector 1, ..., feature vector N], with N total docs.
    """
    point_list = []
    for pid in fv_dict:
        vector_columns = [fv_dict[pid][prop] for prop in fv_dict[pid]]
        point_list.append(vector_columns)
    return np.array(point_list)


if __name__ == "__main__":
    # We need a list of properties for our feature vector.
    feature_list = ['activation_energy', 'volume_per_anion', 
                    'average_bond_length', 'min_cation_distance']
    
    # We also need a database and some documents from it.
    yaml_file = "/Users/dchannah/yaml/paper2.yaml"
    coll = initialize_db(yaml_file)

    # Let's just grab all of the documents in the ionic mobility DB for now.
    docs = coll.find()

    # Now we build our feature vector for all these docs.
    fv = build_feature_vector(feature_list, docs)

    # Now let's do some clustering.
    points = build_point_list(fv)
    points = vq.whiten(points)  # Need to normalize each column first.
    centroids, distortion = vq.kmeans(points, 4)
    indices, _ = vq.vq(points, centroids)

