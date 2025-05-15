import itertools
import random
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
)
from sklearn.model_selection import train_test_split

from network_builder import SimplicialComplexBuilder


def get_edges_from_sets(sets):
    """Return a list of unique edges formed by all 2-element combinations within each input set."""
    edges = set()
    for s in sets:
        for v, u in itertools.combinations(s, 2):
            assert v <= u
            edges.add((v, u))
    return list(edges)


class LinkPredictor:
    MAX_FAILED_ATTEMPTS = 5000

    def __init__(self, sets, test_size=0.3, neg_ratio=1.0, random_state=42):
        """
        sets: list of sets of nodes (each hyperedge)
        test_size: test_size in ML task
        neg_ratio: number of negative samples per positive
        random_state: seed for reproducibility
        """
        self.sets = sets
        self.test_size = test_size
        self.neg_ratio = neg_ratio
        self.random_state = random_state
        self.scaler = StandardScaler()
        self.node_features = []
        self.clique_edges = set()
        self.V = max(max(s) for s in self.sets) + 1
        self.rng = random.Random(self.random_state)

    def sample_negative_edges_sparse_graph(self, num_samples):
        """
        Generate num_samples edges, which are absent
        by selecting random vertices and combining them in edge.
        If threshold reached, then return None, otherwise return list of generated edges.
        """
        neg_edges = set()
        failed_attempts = 0
        while len(neg_edges) < num_samples:
            if failed_attempts == LinkPredictor.MAX_FAILED_ATTEMPTS:
                return None
            v, u = self.rng.sample(range(0, self.V), 2)
            edge = (min(v, u), max(v, u))
            if edge not in self.clique_edges and edge not in neg_edges:
                neg_edges.add(edge)
                failed_attempts = 0
            else:
                failed_attempts += 1
        return list(neg_edges)

    def sample_negative_edges_dense_graph(self, num_samples):
        """
        Generate num_samples edges, which are absent
        by selecting prefix of random permutation of edges.
        """
        edges = list(itertools.combinations(range(self.V), 2))
        self.rng.shuffle(edges)
        neg_edges = set()
        for edge in edges:
            if len(neg_edges) == num_samples:
                break
            if edge not in self.clique_edges and edge not in neg_edges:
                neg_edges.add(edge)
        return list(neg_edges)

    def sample_negative_edges(self, num_samples):
        """
        Generate num_samples edges, which are absent.
        """
        result = self.sample_negative_edges_sparse_graph(num_samples)
        if result is None:
            result = self.sample_negative_edges_dense_graph(num_samples)
        if len(result) < num_samples:
            print("Warning! Generated less negative samples then should be.")
        return result

    def calculate_node_features(self):
        clos1 = np.array(
            [value for _, value in self.network.ClosenessAll(0, 1, False)]
        ).reshape(-1, 1)
        clos2 = np.array(
            [value for _, value in self.network.ClosenessAll(0, 2, False)]
        ).reshape(-1, 1)
        centr_eig1 = np.array(
            [value for _, value in self.network.EigenCentrality(0, 1, False)]
        ).reshape(-1, 1)
        centr_eig2 = np.array(
            [value for _, value in self.network.EigenCentrality(0, 2, False)]
        ).reshape(-1, 1)
        centr_sub1 = np.array(
            [value for _, value in self.network.SubgraphCentrality(0, 2, False)]
        ).reshape(-1, 1)
        centr_sub2 = np.array(
            [value for _, value in self.network.SubgraphCentrality(0, 1, False)]
        ).reshape(-1, 1)
        deg1 = np.array(self.network.DegreeAll(0, 1, False)).reshape(-1, 1)
        deg2 = np.array(self.network.DegreeAll(0, 2, False)).reshape(-1, 1)
        # betw1 = np.array(
        #     [value for _, value in self.network.BetweennessAll(0, 1, False)]
        # ).reshape(-1, 1)
        # betw2 = np.array(
        #     [value for _, value in self.network.BetweennessAll(0, 2, False)]
        # ).reshape(-1, 1)

        # Centr1
        self.node_features = np.hstack((clos1, centr_eig1, centr_sub1, deg1))
        # Centr2
        self.node_features = np.hstack(
            (self.node_features, clos2, centr_eig2, centr_sub2, deg2)
        )

        # Without centr
        # self.node_features = np.zeros((0, 0))

        print(f"Using {self.node_features.shape[1]} node features")

    def extract_features_on_edges(self, edges):
        """
        Return features for every edge
        """
        edge_features = []
        for v, u in edges:
            edge_feature = []
            for i in range(self.node_features.shape[1]):
                edge_feature.append(self.node_features[v][i])
                edge_feature.append(self.node_features[u][i])
            # CN1
            edge_feature.append(self.network.CommonNeighbors([v], [u], 1))
            # CN2
            edge_feature.append(self.network.CommonNeighbors([v], [u], 2))
            # CN3
            # edge_feature.append(self.network.CommonNeighbors([v], [u], 3))

            edge_features.append(np.array(edge_feature))
        return np.vstack(edge_features)

    def make_sets_connected_component(self, sets_train, sets_test):
        """
        By given division sets on train and test return new division,
        in which train sets are augmented with tests set, so it
        will contain all vertices from possible vertices
        """
        missing = set(range(0, self.V))
        for s in sets_train:
            missing -= set(s)
        sets_train_final = sets_train.copy()
        sets_test_final = []
        for s in sets_test:
            if missing & set(s):
                sets_train_final.append(s)
                missing -= set(s)
            else:
                sets_test_final.append(s)

        return sets_train_final, sets_test_final

    def train(self):
        """
        Train logistic regression model on sets features.
        Returns AP and ROC scores.
        """
        # Split sets on train/test
        sets_train, sets_test = train_test_split(
            self.sets, test_size=self.test_size, random_state=self.random_state
        )

        # Add missing sets in train, to make network contain all vertices from 0 to V
        sets_train, sets_test = self.make_sets_connected_component(
            sets_train, sets_test
        )
        print("Divided sets on train/test!")

        # Build network + features on nodes
        self.network = SimplicialComplexBuilder().create_network(sets_train)
        self.calculate_node_features()
        print("Calculated node features!")
        # Get all edges in sets
        self.clique_edges = get_edges_from_sets(self.sets)

        # Get train/test positive edges
        pos_edges_train = get_edges_from_sets(sets_train)
        pos_edges_test = get_edges_from_sets(sets_test)
        # Generate negative samples
        neg_edges = self.sample_negative_edges(
            (len(pos_edges_train) + len(pos_edges_test)) * self.neg_ratio
        )
        # Split negative edges
        neg_edges_train, neg_edges_test = train_test_split(
            neg_edges,
            test_size=len(pos_edges_test)
            / (len(pos_edges_train) + len(pos_edges_test)),
            random_state=self.random_state,
        )
        # Unite pos/neg edges
        edges_train = pos_edges_train + neg_edges_train
        edges_test = pos_edges_test + neg_edges_test
        y_train = np.array([1] * len(pos_edges_train) + [0] * len(neg_edges_train))
        y_test = np.array([1] * len(pos_edges_test) + [0] * len(neg_edges_test))

        # Extract Features
        X_train = self.extract_features_on_edges(edges_train)
        X_test = self.extract_features_on_edges(edges_test)

        # Scale
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        self.model = LogisticRegression(solver="lbfgs", random_state=self.random_state)
        self.model.fit(X_train_scaled, y_train)

        y_pred = self.model.predict_proba(X_test_scaled)[:, 1]

        auc_ap, auc_roc = self.calculate_score(y_test, y_pred)
        print(f"Test AUC-AP: {auc_ap}")
        print(f"Test AUC-ROC: {auc_roc}")
        return (auc_ap, auc_roc)

    def calculate_score(self, y_test, y_prob):
        ap = average_precision_score(y_test, y_prob)
        auc_roc = roc_auc_score(y_test, y_prob)
        return (ap, auc_roc)
