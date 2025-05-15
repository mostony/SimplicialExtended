import numpy as np


class Aggregator:
    """
    Class for dealing with features on higher-dimension structure which
    needs to be propagated down.
    In basic case features on edges need to be push down on vertices.
    """

    def __init__(self, V, aggregate_func=lambda a, b: a + b, default_value=0):
        self.V = V
        self.aggregate_func = aggregate_func
        self.default_value = default_value

    def push_down(self, edge_features, edges):
        assert len(edges) == edge_features.shape[0]
        result = np.full((self.V, edge_features.shape[1]), float(self.default_value))
        for col in range(edge_features.shape[1]):
            for row in range(edge_features.shape[0]):
                feat_value = edge_features[row][col]
                for vertice in edges[row]:
                    assert 0 <= vertice and vertice < self.V
                    result[vertice][col] = self.aggregate_func(
                        result[vertice][col], feat_value
                    )
        return result
