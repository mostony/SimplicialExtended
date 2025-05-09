from typing import List, Tuple


class DatasetLoader:
    def load_dataset(name: str) -> Tuple[List[List[int]], List[int]]:
        """
        Get (hyperedges, labels) by name of dataset
        """
        sets = FileDatasetLoader.load_hyper_edges(
            f"./data/{name}/hyperedges-{name}.txt"
        )
        true_labels = FileDatasetLoader.load_node_labels(
            f"./data/{name}/node-labels-{name}.txt"
        )
        return (sets, true_labels)


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
                true_labels.append(int(line))
        return true_labels


class PytorchDatasetLoader:
    def load_hyper_edges(dataset_name: str) -> List[List[int]]:
        pass

    def load_node_labels(path: str) -> List[int]:
        pass
