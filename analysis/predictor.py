import itertools
import random
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
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
    def __init__(self, sets, test_size=0.3, neg_ratio=1.0, random_state=42):
        """
        hyperedges: list of iterables of nodes (each hyperedge)
        neg_ratio: number of negative samples per positive
        random_state: integer seed for reproducibility
        """
        self.sets = sets
        self.test_size = 0.3
        self.neg_ratio = neg_ratio
        self.random_state = random_state
        self.scaler = StandardScaler()
        self.node_features = []
        self.all_edges = set()
        self.V = max(max(s) for s in self.sets) + 1

    def extract_features_on_edges(self, edges):
        """
        Return features for every edge
        """
        features = []
        for v, u in edges:
            edge_feature = []
            for i in range(self.node_features.shape[1]):
                edge_feature.append(self.node_features[v][i])
                edge_feature.append(self.node_features[u][i])
            
            # CN
            edge_feature.append(sum(int(a) & int(b) for a, b in zip(self.adj1[u], self.adj1[v])))
            # edge_feature.append(sum(int(a) & int(b) for a, b in zip(self.adj2[u], self.adj2[v])))

            # TODO: Bad!!!
            # edge_feature = [(sum(int(a) & int(b) for a, b in zip(self.adj1[u], self.adj1[v])))]
            features.append(np.array(edge_feature))
        return np.vstack(features)

    def sample_negative_edges(self, num_samples):
        random.seed(self.random_state)
        neg_edges = set()
        iter = 0
        while len(neg_edges) < num_samples:
            v, u = random.sample(range(0, self.V), 2)
            assert u != v
            if v > u:
                v, u = u, v
            if (v, u) not in self.all_edges:
                neg_edges.add((v, u))
                iter = 0
            if iter == 100:
                raise RuntimeError("Can't find negative edge")
        return list(neg_edges)

    def sample_negative_hyper_edges(self, num_samples):
        random.seed(self.random_state)
        neg_edges = set()
        iter = 0
        while len(neg_edges) < num_samples:
            v, u = random.sample(range(0, self.V), 2)
            assert u != v
            if v > u:
                v, u = u, v
            if (v, u) not in self.all_edges:
                neg_edges.add((v, u))
                iter = 0
            if iter == 100:
                raise RuntimeError("Can't find negative edge")
        return list(neg_edges)

    def calculate_node_features(self):
        # clos1 = np.array(
        #     [value for _, value in self.network.ClosenessAll(0, 1, False)]
        # ).reshape(-1, 1)
        # clos2 = np.array(
        #     [value for _, value in self.network.ClosenessAll(0, 2, False)]
        # ).reshape(-1, 1)
        # centr_eig1 = np.array(
        #     [value for _, value in self.network.EigenCentrality(0, 1, False)]
        # ).reshape(-1, 1)
        # centr_eig2 = np.array(
        #     [value for _, value in self.network.EigenCentrality(0, 2, False)]
        # ).reshape(-1, 1)
        # centr_sub1 = np.array(
        #     [value for _, value in self.network.SubgraphCentrality(0, 2, False)]
        # ).reshape(-1, 1)
        # centr_sub2 = np.array(
        #     [value for _, value in self.network.SubgraphCentrality(0, 1, False)]
        # ).reshape(-1, 1)
        # deg1 = np.array(self.network.DegreeAll(0, 1, True)).reshape(-1, 1)
        # deg2 = np.array(self.network.DegreeAll(0, 2, True)).reshape(-1, 1)
        betw1 = np.array(
            [value for _, value in self.network.BetweennessAll(0, 1, False)]
        ).reshape(-1, 1)
        betw2 = np.array(
            [value for _, value in self.network.BetweennessAll(0, 2, False)]
        ).reshape(-1, 1)

        # self.node_features = np.hstack((clos1, centr_eig1, centr_sub1))
        # self.node_features = np.hstack((self.node_features, clos2, centr_eig2, centr_sub2))
        # self.node_features = np.hstack((self.node_features, deg1, deg2, betw1, betw2))
        # self.node_features = betw1
        self.node_features = np.zeros((0, 0)) 
        self.adj1 = self.network.AdjacencyMatrix(0, -1, 1, False)
        # self.adj2 = self.network.AdjacencyMatrix(0, -1, 2, False)
        print(f"Using {self.node_features.shape[1]} node features")
        return
    
    def make_sets_connected_component(self, sets_train, sets_test):
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
        """Train logistic regression model on sets features
        Returns AP score.
        """
        # Split sets on train/test
        sets_train, sets_test = train_test_split(
            self.sets, test_size=self.test_size, random_state=self.random_state
        )
        # Add missing sets, to make network contain all vertices from 0 to V
        sets_train, sets_test = self.make_sets_connected_component(sets_train, sets_test)
        
        # Build network + features on nodes
        self.network = SimplicialComplexBuilder().create_network(sets_train)
        self.calculate_node_features()
        print("Calculated node features!")
        # Get all edges in sets
        self.all_edges = get_edges_from_sets(self.sets)
        # Generate negative samples
        neg_edges = self.sample_negative_edges(len(self.all_edges) * self.neg_ratio)
        # Split negative edges
        neg_edges_train, neg_edges_test = train_test_split(
            neg_edges, test_size=self.test_size, random_state=self.random_state
        )
        # Get train/test positive edges
        pos_edges_train = get_edges_from_sets(sets_train)
        pos_edges_test = get_edges_from_sets(sets_test)
        # Unite pos/neg edges
        edges_train = pos_edges_train + neg_edges_train
        edges_test = pos_edges_test + neg_edges_test
        y_train = np.array([1] * len(pos_edges_train) + [0] * len(neg_edges_train))
        y_test = np.array([1] * len(pos_edges_test) + [0] * len(neg_edges_test))
        # Get Features
        X_train = self.extract_features_on_edges(edges_train)
        X_test = self.extract_features_on_edges(edges_test)

        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        print("Prepared edge features")
        self.model = LogisticRegression(
            solver="lbfgs", random_state=self.random_state
        )
        self.model.fit(X_train_scaled, y_train)

        train_probs = self.model.predict_proba(X_train_scaled)[:, 1]
        test_probs = self.model.predict_proba(X_test_scaled)[:, 1]
        yes = np.ones(y_test.shape)
        print(f"Train AUC-ROC: {roc_auc_score(y_train, train_probs):.3f}")
        print(f"Train AP: {average_precision_score(y_train, train_probs):.3f}")
        print(f"Test AUC-ROC: {roc_auc_score(y_test, test_probs):.3f}")
        print(f"Test AP: {average_precision_score(y_test, test_probs):.3f}")

        print(f"Yes Test AUC-ROC: {roc_auc_score(y_test, yes):.3f}")
        print(f"Yes Test AP: {average_precision_score(y_test, yes):.3f}")
        score = average_precision_score(y_test, test_probs)
        return score

    def predict_proba(self, hyperedge):
        """Return probability that hyperedge will appear in future"""
        h_set = set(hyperedge)
        feats = np.array(self._extract_features(h_set)).reshape(1, -1)
        feats_scaled = self.scaler.transform(feats)
        return self.model.predict_proba(feats_scaled)[0, 1]
