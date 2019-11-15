import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris
iris = load_iris()

X = pd.DataFrame(iris.data, columns=iris.feature_names)
y = pd.Categorical.from_codes(iris.target, iris.target_names)

print X.head()


scaler = StandardScaler()
X = scaler.fit_transform(X)

pca = PCA(n_components=4)

principal_components = pca.fit_transform(X)new_X = pd.DataFrame(data = principal_components, columns = ['PC1', 'PC2', 'PC3', 'PC4'])
