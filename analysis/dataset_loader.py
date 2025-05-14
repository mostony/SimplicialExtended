from typing import List, Tuple


class DatasetLoader:

    def load_dataset(name: str, load_labels=True):
        """
        Get (hyperedges, labels) by name of dataset
        """
        sets = FileDatasetLoader.load_hyper_edges(
            f"./data/{name}/hyperedges-{name}.txt"
        )

        # Get all not isolated vertices
        not_isolated = set()
        for s in sets:
            for v in s:
                not_isolated.add(v)

        # Compress vertices in [0, V)
        sorted_values = sorted(not_isolated)
        value_to_compressed = {val: idx for idx, val in enumerate(sorted_values)}

        filtered_sets = []
        for row in sets:
            s = row
            for i in range(len(s)):
                s[i] = value_to_compressed[s[i]]

            if len(s) == 1:
                s.append(s[0])
            assert len(s) > 0
            filtered_sets.append(s)

        if load_labels:
            true_labels = FileDatasetLoader.load_node_labels(
                f"./data/{name}/node-labels-{name}.txt"
            )
            filtered_labels = [true_labels[i] for i in sorted(not_isolated)]
            return filtered_sets, filtered_labels

        return filtered_sets

    @staticmethod
    def _add_loops(sets: List[List[int]], V: int) -> List[List[int]]:
        """
        Adds for every vertex loop to deal with isolated vertices
        """
        for v in range(V):
            sets.append([v, v])
        return sets


class FileDatasetLoader:
    MAX_EDGE_SIZE = 10

    def load_hyper_edges(path: str) -> List[List[int]]:
        """
        File by given path contains set of hyperedges. Each hyperedge on distinct line, vertices separated by comma.
        Returns list of hyper_edges. Vertices in each hyper_edge are sorted and numeration from zero.
        Hyperedges with size >= MAX_EDGE_SIZE are ignored.
        """
        hyper_edges = []

        with open(path, "r") as file:
            for line in file:
                edge = list(map(int, line[:-1].split(",")))
                edge.sort()
                edge = list(map(lambda x: x - 1, edge))
                assert min(edge) >= 0 and len(edge) == len(set(edge))
                if len(edge) <= FileDatasetLoader.MAX_EDGE_SIZE:
                    hyper_edges.append(edge)

        return hyper_edges

    def load_node_labels(path: str) -> List[int]:
        """
        File by given path contains labels for each vertices. Each label on distinct line.
        Returns true label for every node in dataset.
        """
        true_labels = []
        with open(path, "r") as file:
            for line in file:
                true_labels.append(int(line) - 1)
        return true_labels


class PytorchDatasetLoader:
    def load_hyper_edges(dataset_name: str) -> List[List[int]]:
        pass

    def load_node_labels(path: str) -> List[int]:
        pass
