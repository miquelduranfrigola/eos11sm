import os
import csv
import joblib
import numpy as np
from FPSim2 import FPSim2Engine
from sklearn.linear_model import LogisticRegressionCV
from sklearn.utils import shuffle

root = os.path.abspath(os.path.dirname(__file__))
checkpoints_dir = os.path.join(root, '..', '..', '..', 'checkpoints')

FP_NAMES = [
    "morgan_2048",
    "maccs_166",
    "rdkit_2048",
    "atompair_2048",
    "pattern_2048"
]


def get_X(smiles_list, skip_identity=True):
    R = []
    headers = []
    for fp_name in FP_NAMES:
        print(f"Processing fingerprint: {fp_name}")
        db_path = os.path.join(
            checkpoints_dir,
            f"{fp_name}.h5"
        )
        fpe = FPSim2Engine(db_path)
        k1s = []
        k3s = []
        if skip_identity:
            k = 4
            idx_1 = 1
            idx_3 = 3
        else:
            k = 3
            idx_1 = 0
            idx_3 = 2
        for smiles in smiles_list:
            result = fpe.top_k(smiles, k=k, threshold=0.0, metric='tanimoto', n_workers=1)
            k1s += [float(result[idx_1][1])]
            k3s += [float(result[idx_3][1])]
        R += [k1s, k3s]
        headers += [f"{fp_name}_k1", f"{fp_name}_k3"]
    X = np.array(R).T
    return X
        

def get_positive_X_y():
    with open(os.path.join(root, "..", "results", "00_abx_compounds.csv"), "r") as f:
        reader = csv.reader(f)
        next(reader)
        smiles_list = []
        for r in reader:
            smiles_list.append(r[0])
    print(smiles_list[:5])
    y = np.array([1]*len(smiles_list))
    X = get_X(smiles_list, skip_identity=True)
    return X, y


def get_negative_X_y():
    with open(os.path.join(root, "..", "results", "03_non_abx_compounds.csv"), "r") as f:
        reader = csv.reader(f)
        next(reader)
        smiles_list = []
        for r in reader:
            smiles_list.append(r[0])
    print(smiles_list[:5])
    y = np.array([0]*len(smiles_list))
    X = get_X(smiles_list, skip_identity=False)
    return X, y


X1, y1 = get_positive_X_y()
X0, y0 = get_negative_X_y()

X = np.vstack([X1, X0])
y = np.hstack([y1, y0])
X, y = shuffle(X, y, random_state=42)

print(X, y)

model = LogisticRegressionCV(class_weight='balanced', max_iter=1000, scoring='roc_auc', cv=5)
model.fit(X, y)
print(model.scores_[1])


joblib.dump(model, os.path.join(root, "..", "..", "..", "checkpoints", "logistic_regression.joblib"))
