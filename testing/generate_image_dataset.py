from collections import Counter
from imblearn.datasets import fetch_datasets
from sklearn.neighbors import kneighbors_graph
import numpy as np
import time

import os


def FetchImageDataset():
    print(f"Fetching dataset...")
    result = fetch_datasets()["letter_img"].data
    print(f"Finished fetching dataset!")
    return result


def GenerateImageDataset(n, k, random_seed=228):
    """
    Generate subsample with size = n of letter_img dataset.
    Then build kneighbors_graph and symmetrize it

    Args:
        n : size of subsample
        k : number of nearest neighbours for every point
        random_seed: random seed to choose subsample
    """
    if not hasattr(GenerateImageDataset, "dataset"):
        GenerateImageDataset.dataset = FetchImageDataset()

    letter_img_dataset = GenerateImageDataset.dataset
    N = letter_img_dataset.shape[0]

    # TODO: add some condition about size of subsample

    subsample = letter_img_dataset[np.random.choice(N, n, replace=False), :]

    print(f"Generating graph for n={n}, k={k}...")
    start_time = time.time()

    A = kneighbors_graph(subsample, n_neighbors=k, n_jobs=-1)
    A = ((A + A.T) > 0).astype(int).A

    end_time = time.time()

    with open(f"testing/times{k}.txt", "a+") as f:
        f.write(str(end_time - start_time) + " ")

    print(f"Finished generating graph for n={n}, k={k}!")
    print(f"Elapsed time of generating graph : {end_time - start_time} s")

    print(f"Saving file...")

    np.savetxt(f"testing/images.txt", A, fmt="%i", delimiter=" ")
    print(f"Finished saving file!")

    # return A


import sys

if __name__ == "__main__":
    assert len(sys.argv) == 4
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    random_seed = int(sys.argv[3])
    GenerateImageDataset(n, k, random_seed)
