# imports
import os
import csv
import sys
import joblib
import numpy as np
from FPSim2 import FPSim2Engine

# variables
FP_NAMES = [
    "morgan_2048",
    "maccs_166",
    "rdkit_2048",
    "atompair_2048",
    "pattern_2048"
]

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))
checkpoints_dir = os.path.join(root, '..', '..', 'checkpoints')

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)
    smiles_list = [r[0] for r in reader]

lr_model = joblib.load(os.path.join(checkpoints_dir, "logistic_regression.joblib"))

def get_X(smiles_list):
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

# run
X = get_X(smiles_list)
outputs = lr_model.predict_proba(X)[:, 1].tolist()

# check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["abx_score"])
    for o in outputs:
        writer.writerow([o])
