from pylab import plot, show, scatter
from analyze_db import build_db_query, initialize_db, build_feature_vector
from matplotlib import cm
from sklearn import cluster, linear_model, cross_validation
import scipy.cluster.vq as vq
import numpy as np


def regress_feature(regression_model=linear, cluster, metric, obs):
    """regress_feature
    Perform regression of values of some metric against a observed data.
    :param regression_model: Regression model to use (default: linear)
    :param metric: A list of metric measurements for each point in a cluster
    :param obs: A list of observations for each point in the cluster
    :return: Mean absolute error of regression (metric vs. obs)
    """
    if len(metric) != len(obs):
        print "Observations and feature vector item length are not the same."
        return None
    else:
        cv = cross_validation.ShuffleSplit(len(obs), n_iter=10, test_size=0.1,
                                                   random_state=0)
        scores = cross_validation.cross_val_score(linear, metric, obs, 
                                                  scoring='mean_absolute_error')
        return np.abs(np.mean(scores))


def build_traj_pt_list(fv_dict):
    """build_traj_pt_list
    Builds a point list where each point is a set of props for an NEB image.
    :param fv_dict: Dictionary of traj-prop feature vectors.
    :return: [[trp_A^img1, ..., trp_n^img1], ... [trp_A^imgN..., trp_n^imgN]]
    """
    # Feature vector dict items are a list of lists now.
    point_list = []
    for neb_traj in fv_dict:
        traj_array = []
        for prop in fv_dict[neb_traj]:
            traj_array.append(fv_dict[neb_traj][prop])
        for i in range(len(traj_array[0])):
            point_list.append(np.array(traj_array)[:,i])
    return np.array(point_list)


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
    feature_list = ['min_anion_distance', 'volume_per_anion', 
                    'average_bond_length', 'min_cation_distance']

    # Also testing a trajectory-resolved feature vector
    traj_fl = ['effective_coordination', 'min_cat_dist',
               'min_anion_dist']
    
    # We also need a database and some documents from it.
    yaml_file = "/Users/dchannah/yaml/paper2.yaml"
    coll = initialize_db(yaml_file)

    # docs = coll.find({"cation_type": "Mg"})  # Mg-ion cathodes only
    docs = coll.find()  # All documents
    
    # Now we build our feature vector for all these docs.
    fv = build_feature_vector(feature_list, docs)
    # fv_fl = build_feature_vector(traj_fl, docs)

    # From our DB feature vector, need to put data in array form
    # points = build_traj_pt_list(fv_fl)
    points = build_point_list(fv)
    points = vq.whiten(points)

    # Now we get the k means using KMeans++
    k = 2
    kmeans = cluster.KMeans(n_clusters=k)
    kmeans.fit(points)

    labels, centroids = kmeans.labels_, kmeans.cluster_centers_

    # Now we plot
    cmap = cm.jet
    for i in range(k):
        # plot points
        points_in_cluster = points[np.where(labels==i)]
        scatter(points_in_cluster[:,0], points_in_cluster[:,1], color=cmap(i*50))

        # plot centroids
        scatter(centroids[i,0], centroids[i,1], s=80, marker='x')
    show()

