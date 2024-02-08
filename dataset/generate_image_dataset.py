from collections import Counter
from imblearn.datasets import fetch_datasets
from sklearn.neighbors import kneighbors_graph
import numpy as np

letter_img_dataset = fetch_datasets()["letter_img"]

# print(letter_img_dataset.data[:10])
# print(letter_img_dataset.target[:10])
# print(sorted(Counter(letter_img_dataset.target).items()))

# it takes around 20 min
for n_neighbors in [5, 10, 15, 20]:
    print(f"Generating graph for n_neighbors = {n_neighbors}...")
    A = kneighbors_graph(letter_img_dataset.data, n_neighbors=n_neighbors, n_jobs=-1)
    A = ((A + A.T) > 0).astype(int).A
    np.savetxt(f"images{n_neighbors}.txt", A, fmt="%i", delimiter=" ")
