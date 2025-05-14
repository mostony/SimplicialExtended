import numpy as np


class Aggregator:
    # TODO: make default value 
    def __init__(self, V, aggregate_func=sum):
        self.V = V
        self.aggregate_func = aggregate_func

    def push_down(self, edge_feature, edges):
        assert len(edges) == len(edge_feature)
        result = np.zeros((self.V, 1))
        for i, feature in enumerate(edge_feature):
            for vertice in edges[i]:
                assert vertice < self.V
                result[vertice] = self.aggregate_func(result[vertice], feature)
        return result
