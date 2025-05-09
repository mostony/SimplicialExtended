from simpl_loader import load_simpl_module
from itertools import chain, combinations
from typing import Any, List
from abc import ABC, abstractmethod

simpl = load_simpl_module()


# TODO: cache?
class Builder(ABC):
    def __init__(self, sets: List[List[int]], weighted=True):
        self.sets = sets
        self.edges = {}
        self.structure = None
        self.weighted = weighted

    def build(self):
        self.collect_edges()
        self.create_structure()
        if self.weighted:
            self.update_weights()
        return self.structure

    @abstractmethod
    def collect_edges(self):
        pass

    @abstractmethod
    def create_structure(self):
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

    def create_structure(self):
        self.structure = simpl.SimplicialComplex()
        for edge in self.sets:
            self.structure.AddSimplex(edge)

    def update_weights(self):
        for subset, weight in self.edges.items():
            self.structure.UpdateWeight(list(subset), weight)

    def _get_all_subsets(self, s: List[int]) -> List[List[int]]:
        return list(chain.from_iterable(combinations(s, r) for r in range(1, len(s) + 1))) 

class HyperGraphBuilder(Builder):
    def collect_edges(self):
        for edge in self.sets:
            edge = tuple(edge)
            self.edges[edge] = self.edges.get(edge, 0) + 1

    def create_structure(self):
        self.structure = simpl.HyperGraph()
        for edge in self.sets:
            self.structure.AddEdge(edge)

    def update_weights(self):
        for edge, weight in self.edges.items():
            self.structure.UpdateWeight(list(edge), weight)


class GraphBuilder(Builder):
    def collect_edges(self):
        for hyper_edge in self.sets:
            for v in hyper_edge:
                for u in hyper_edge:
                    if v < u:
                        edge = (v, u)
                        self.edges[edge] = self.edges.get(edge, 0) + 1

    def create_structure(self):
        self.structure = simpl.Graph()
        for (v, u) in self.edges:
            self.structure.AddEdge(v, u)

    def update_weights(self):
        for (v, u), weight in self.edges.items():
            self.structure.UpdateWeight((v, u), weight)