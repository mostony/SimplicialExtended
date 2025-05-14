from simpl_loader import load_simpl_module
from itertools import chain, combinations
from typing import Any, List
from abc import ABC, abstractmethod

simpl = load_simpl_module()


# TODO: cache?
class Builder(ABC):

    def create_network(self, sets: List[List[int]], weighted=True):
        self.edges = {}
        self.vertices = {}
        self.network = None
        self.weighted = weighted
        self.sets = sets

        self.collect_edges()
        self.init_network()
        if self.weighted:
            self.update_weights()
        return self.network

    @abstractmethod
    def collect_edges(self):
        pass

    def collect_vertices(self):
        for edge, w in self.edges.items():
            for vertice in edge:
                self.vertices[vertice] = self.vertices.get(vertice, 0) + w

    @abstractmethod
    def init_network(self):
        pass

    @abstractmethod
    def update_weights(self):
        pass


class SimplicialComplexBuilder(Builder):
    def collect_edges(self):
        for edge in self.sets:
            for subset in self._get_all_subsets(edge):
                subset = tuple(subset)
                self.edges[subset] = self.edges.get(subset, 0) + 1

    def init_network(self):
        self.network = simpl.SimplicialComplex()
        for edge in self.sets:
            self.network.AddSimplex(edge)

    def update_weights(self):
        for subset, weight in self.edges.items():
            self.network.UpdateWeight(list(subset), weight)

    def _get_all_subsets(self, s: List[int]) -> List[List[int]]:
        return list(
            chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))
        )


class HyperGraphBuilder(Builder):
    def collect_edges(self):
        for edge in self.sets:
            edge = tuple(edge)
            self.edges[edge] = self.edges.get(edge, 0) + 1

    def init_network(self):
        self.network = simpl.HyperGraph()
        for edge in self.sets:
            self.network.AddEdge(edge)

    def update_weights(self):
        # Update edges
        for edge, weight in self.edges.items():
            self.network.UpdateWeight(list(edge), weight)
        # Update vertices
        self.collect_vertices()
        for vertice, weight in self.vertices.items():
            self.network.UpdateWeight([vertice], weight)


class GraphBuilder(Builder):
    def collect_edges(self):
        for hyper_edge in self.sets:
            for pos_v, v in enumerate(hyper_edge):
                for pos_u, u in enumerate(hyper_edge):
                    if pos_v < pos_u:
                        edge = (v, u)
                        self.edges[edge] = self.edges.get(edge, 0) + 1

    def init_network(self):
        self.network = simpl.Graph()
        for v, u in self.edges:
            self.network.AddEdge(v, u)

    def update_weights(self):
        for (v, u), weight in self.edges.items():
            self.network.UpdateWeight((v, u), weight)

        # Update vertices
        self.collect_vertices()
        for vertice, weight in self.vertices.items():
            self.network.UpdateWeight([vertice], weight)
