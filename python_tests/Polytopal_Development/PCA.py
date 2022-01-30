import numpy as np
import matplotlib.pyplot as plt


from sklearn.decomposition import PCA

# creating the dataset
x=np.array([[1,3,0],[2,2,0],[3,4,0],[4,2,0],[5,9,0]])
plt.scatter(x[:,0], x[:,1])
print(x)
pca = PCA(n_components=3).fit(x)
print(pca.components_)
pca = PCA(n_components=2).fit(x)
pca_x = pca.transform(x)
print(pca.components_)
x2 = pca.inverse_transform(pca_x)
plt.scatter(x[:,0], x[:,1], alpha=0.3)
plt.scatter(x2[:,0], x2[:,1])
plt.show()
print(pca_x)
