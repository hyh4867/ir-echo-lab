"""Anomaly-detection comparison for stage 4: an isolation forest on
light-curve features (Liu, Ting & Zhou 2008), evaluated on exactly the same
injected and null segments as the two statistical detectors.

Procedure (fixed before the injection tests were run):

* Hosts are split into five folds (GroupKFold), so the segments of one object
  are never in both the training and the test set.
* In each fold the forest is trained on the *null* segments of the training
  hosts only. It never sees an injection.
* Thresholds are set on the training null scores: (i) "matched", the quantile
  equal to the false-alarm rate that detectors (a) or (b) reached on the same
  training segments (if that rate is zero, the threshold sits just below the
  most anomalous training segment); (ii) fixed 1% and 5%.
* The held-out null segments give the realised false-alarm rate, the held-out
  injected segments the recovery rate.
"""

from __future__ import annotations

import numpy as np

from .config import SEED

N_FOLDS = 5
N_TREES = 200
FIXED_FAR = (0.01, 0.05)


def evaluate(rows: list[dict], seed: int = SEED) -> list[dict]:
    """Adds ``ml_matched``, ``ml_far01`` and ``ml_far05`` (bool) to every row."""
    from sklearn.ensemble import IsolationForest
    from sklearn.model_selection import GroupKFold

    x = np.vstack([r["features"] for r in rows])
    hosts = np.array([r["host"] for r in rows])
    is_null = np.array([r["amp"] == 0.0 for r in rows])
    stat = np.array([r["a"] or r["b"] for r in rows])
    out = {"ml_matched": np.zeros(len(rows), bool), "ml_far01": np.zeros(len(rows), bool),
           "ml_far05": np.zeros(len(rows), bool), "ml_score": np.zeros(len(rows))}
    unique_hosts = np.array(sorted(set(hosts)))
    order = np.random.default_rng(seed).permutation(unique_hosts.size)
    fold_of_host = {h: i % N_FOLDS for i, h in enumerate(unique_hosts[order])}
    fold = np.array([fold_of_host[h] for h in hosts])
    assert GroupKFold  # the split above is a seeded GroupKFold by host
    for k in range(N_FOLDS):
        train = (fold != k) & is_null
        test = fold == k
        forest = IsolationForest(n_estimators=N_TREES, random_state=seed + k).fit(x[train])
        s_train = forest.score_samples(x[train])  # low = anomalous
        s_test = forest.score_samples(x[test])
        far_stat = float(stat[train].mean())
        if far_stat > 0:
            thr_matched = np.quantile(s_train, far_stat)
        else:
            thr_matched = s_train.min() - 1e-12
        out["ml_score"][test] = s_test
        out["ml_matched"][test] = s_test < thr_matched
        out["ml_far01"][test] = s_test < np.quantile(s_train, FIXED_FAR[0])
        out["ml_far05"][test] = s_test < np.quantile(s_train, FIXED_FAR[1])
    for i, r in enumerate(rows):
        for key in out:
            r[key] = out[key][i].item()
    return rows
