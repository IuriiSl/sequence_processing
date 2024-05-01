import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from multiprocessing import Pool


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(self, n_estimators=10, max_depth=None, max_features=None, random_state=42):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit_loop_func(self, args):
        X, y, i = args
        np.random.seed(self.random_state + i)
        feat_idx = np.random.choice(X.shape[1], size=self.max_features, replace=False).tolist()

        sample_indices = np.random.choice(range(X.shape[0]), size=X.shape[0], replace=True)
        X_sampled = X[sample_indices][:, feat_idx]
        y_sampled = y[sample_indices]

        clf = DecisionTreeClassifier(max_depth=self.max_depth, max_features=self.max_features,
                                     random_state=self.random_state)
        clf.fit(X_sampled, y_sampled)
        return clf, feat_idx

    def fit(self, X, y, n_jobs: int = 1):
        self.classes_ = sorted(np.unique(y))

        with Pool(n_jobs) as pool:
            self.trees, self.feat_ids_by_tree = zip(*pool.map(self.fit_loop_func,
                                                              [(X, y, i) for i in range(self.n_estimators)]))

        return self

    def predict_proba_loop_func(self, args):
        X, i, tree = args
        X_sample = X[:, self.feat_ids_by_tree[i]]
        return tree.predict_proba(X_sample)

    def predict_proba(self, X, n_jobs: int = 1):
        with Pool(n_jobs) as pool:
            probability = sum(pool.map(self.predict_proba_loop_func,
                                       [(X, i, tree) for i, tree in enumerate(self.trees)]))

        mean_proba = probability / len(self.trees)
        return mean_proba

    def predict(self, X, n_jobs: int = 1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
