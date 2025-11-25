import os
import csv
from tqdm import tqdm
import numpy as np
import pandas as pd

root = os.path.join(os.path.dirname(__file__))

POS_NEG_RATIO = 100

da = pd.read_csv(os.path.join(root, "..", "results", "02_abx_properties.csv"))
db = pd.read_csv(os.path.join(root, "..", "results", "02_background_properties.csv"))

X = np.zeros((da.shape[0], db.shape[0]), dtype=int)

for i, r0 in tqdm(enumerate(da.iterrows())):
    r = r0[1]
    min_mw = r["mw"] - 50
    max_mw = r["mw"] + 50
    min_logp = r["logp"] - 0.5
    max_logp = r["logp"] + 0.5
    min_nrot = r["nrot"] - 1
    max_nrot = r["nrot"] + 1
    min_nhacc = r["nhacc"] - 1
    max_nhacc = r["nhacc"] + 1
    min_nhdon = r["nhdon"] - 1
    max_nhdon = r["nhdon"] + 1
    for j, r1 in enumerate(db.iterrows()):
        r_bg = r1[1]
        c = 0
        if min_mw <= r_bg["mw"] <= max_mw:
            c += 1
        if min_logp <= r_bg["logp"] <= max_logp:
            c += 1
        if min_nrot <= r_bg["nrot"] <= max_nrot:
            c += 1
        if min_nhacc <= r_bg["nhacc"] <= max_nhacc:
            c += 1
        if min_nhdon <= r_bg["nhdon"] <= max_nhdon:
            c += 1
        X[i, j] = c


def sample_idxs():
    sampled_idxs = set()
    expected_num = min(da.shape[0] * POS_NEG_RATIO, db.shape[0])
    for c in [5, 4, 3, 2, 1, 0]:
        print(f"Sampling for c={c}, sampled so far: {len(sampled_idxs)} / {expected_num}")
        for _ in range(10):
            for i in range(X.shape[0]):
                idxs = np.where(X[i, :] == c)[0]
                np.random.shuffle(idxs)
                sampled_idxs.update(idxs[:POS_NEG_RATIO].tolist())
                if len(sampled_idxs) >= expected_num:
                    return sorted(sampled_idxs)
    return sorted(sampled_idxs)

sampled_idxs = sample_idxs()

smiles_list = db["smiles"].tolist()
neg_smiles_list = [smiles_list[i] for i in sampled_idxs]
with open(os.path.join(root, "..", "results", "03_non_abx_compounds.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["smiles"])
    for smi in neg_smiles_list:
        writer.writerow([smi])
