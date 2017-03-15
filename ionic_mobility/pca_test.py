import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from analyze_db import initialize_db, build_db_query, build_feature_vector 
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as skpca
from scipy.spatial import ConvexHull

"""
Here we attempt principle component analysis on feature vectors comprised of 
various properties of NEB trajectories.
"""

# We need a list of properties for our feature vector.
feature_list = ['activation_energy', 'volume_per_anion', 'average_bond_length',
                'min_cation_distance', 'cation_type']

# Now, we need to initialize the Mongo collection containing our NEB data.
yaml_file = "/Users/dchannah/yaml/paper2.yaml"
coll = initialize_db(yaml_file)

# Now we get the documents which use multivalent working ions.
q = """{\"$or\": [{\"cation_type\": \"Mg\"}, {\"cation_type\": \"Ca\"}, 
    {\"cation_type\": \"Zn\"}]}"""

search_mv = build_db_query(ucq=q, custom=True)
mv_docs = coll.find(search_mv)

# Now we generate the feature vector for each MV document and make a dataframe.
fv_dict = build_feature_vector(feature_list, mv_docs)
df = pd.DataFrame.from_dict(fv_dict, orient='index')

# Setup x and y data.
no_of_measurements = len(feature_list) - 1  # All measurements + 1 energy obs.
x = df.ix[:,0:no_of_measurements].values
y = df.ix[:,no_of_measurements].values

# Transform the data
sklearn_pca = skpca(n_components=2)  # Two principal components.
x_std = StandardScaler().fit_transform(x)
y_sklearn = sklearn_pca.fit_transform(x_std)
coefficients = sklearn_pca.components_.T

with plt.style.context('seaborn-whitegrid'):
    pca1 = 0
    pca2 = 1
    plt.figure(figsize=(6,6))
    for lab, col in zip(('Ca', 'Mg', 'Zn'), ('blue', 'red', 'green')):
        plt.scatter(y_sklearn[y == lab, pca1], y_sklearn[y == lab, pca2],
                    label=lab, c=col)
        points = np.column_stack((y_sklearn[y == lab, pca1], y_sklearn[y == lab, pca2]))
        hull1 = ConvexHull(points)  # A convex hull helps to visualize spread.
        for simplex in hull1.simplices:
            plt.plot(points[simplex, 0], points[simplex, 1], 'k-', c=col)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.show()
