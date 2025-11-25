import os
import csv
from rdkit import Chem
from rdkit.Chem import Descriptors

root = os.path.join(os.path.dirname(__file__))


def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    properties = {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "nrot": Descriptors.NumRotatableBonds(mol),
        "nhacc": Descriptors.NumHAcceptors(mol),
        "nhdon": Descriptors.NumHDonors(mol),
    }
    return properties

abx_smiles = []
with open(os.path.join(root, "..", "results", "00_abx_compounds.csv"), "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        abx_smiles.append(r[0])

bkg_smiles = []
with open(os.path.join(root, "..", "results", "00_background.csv"), "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        bkg_smiles.append(r[0])

header = ["smiles", "mw", "logp", "nrot", "nhacc", "nhdon"]
R = []
for smi in abx_smiles:
    props = calculate_properties(smi)
    if props:
        R.append([smi, props["mw"], props["logp"], props["nrot"], props["nhacc"], props["nhdon"]])
with open(os.path.join(root, "..", "results", "02_abx_properties.csv"), "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(R)

R = []
for smi in bkg_smiles:
    props = calculate_properties(smi)
    if props:
        R.append([smi, props["mw"], props["logp"], props["nrot"], props["nhacc"], props["nhdon"]])
with open(os.path.join(root, "..", "results", "02_background_properties.csv"), "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(R)
