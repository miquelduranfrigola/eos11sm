import os
import csv
import datamol as dm
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors
from standardiser import standardise

root = os.path.join(os.path.dirname(__file__))

dm.disable_rdkit_log()
DATASET_FILE = os.path.join(root, "..", "data", "20251113_List of antibiotics_2025(Combined).csv")


def preprocess_with_datamol(smiles):
    mol = dm.to_mol(smiles)
    mol = dm.fix_mol(mol)
    mol = dm.sanitize_mol(mol)
    smiles = dm.to_smiles(mol)
    return smiles

def preprocess_with_standardiser(smiles):
    mol = Chem.MolFromSmiles(smiles)
    try:
        mol = standardise.run(mol)
    except:
        mol = None
    if mol is None:
        return None
    smiles = Chem.MolToSmiles(mol)
    return smiles

def preprocess_with_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    smiles = Chem.MolToSmiles(mol)
    return smiles

def get_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    inchikey = Chem.MolToInchiKey(mol)
    return inchikey

def preprocess(smiles):
    smiles_0 = preprocess_with_datamol(smiles)
    if smiles_0 is not None:
        smiles_1 = preprocess_with_standardiser(smiles_0)
        if smiles_1 is not None:
            return smiles_1, get_inchikey(smiles_1)
        else:
            print("Could not process with standardiser:", smiles_0)
            return smiles_0, get_inchikey(smiles_0)
    else:
        print("Could not process with datamol:", smiles)
        smiles_0 = preprocess_with_standardiser(smiles)
        if smiles_0 is not None:
            return smiles_0, get_inchikey(smiles_0)
        else:
            print("Could not process with standardiser either:", smiles)
            smiles_1 = preprocess_with_rdkit(smiles)
            if smiles_1 is not None:
                return smiles_1, get_inchikey(smiles_1)
            else:
                print("Could not process with rdkit either:", smiles)
                return None, None


with open(DATASET_FILE, "r") as f:
    reader = csv.reader(f)
    header = next(reader)
    smiles_list = []
    std_smiles_list = []
    std_inchikey_list = []
    for row in tqdm(reader):
        smiles = row[2]
        std_smiles, std_inchikey = preprocess(smiles)
        if std_smiles is not None:
            mol = Chem.MolFromSmiles(std_smiles)
            mw = Descriptors.MolWt(mol)
            if 100 <= mw <= 1500:
                std_smiles_list.append(std_smiles)
                std_inchikey_list.append(std_inchikey)
            else:
                print("Molecular weight out of range (100-1000):", std_smiles, mw)

ik2smi = {}
for smi, ik in zip(std_smiles_list, std_inchikey_list):
    ik2smi[ik] = smi

with open(os.path.join(root, "..", "results", "00_abx_compounds.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["smiles", "inchikey"])
    for ik, smi in ik2smi.items():
        writer.writerow([smi, ik])


with open(os.path.join(root, "..", "data", "chembl_36_chemreps.txt"), "r") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    neg_smiles_list = []
    neg_inchikey_list = []
    for i, row in tqdm(enumerate(reader)):
        smiles = row[1]
        std_smiles, std_inchikey = preprocess(smiles)
        if std_inchikey in ik2smi:
            continue
        if std_smiles is not None:
            mol = Chem.MolFromSmiles(std_smiles)
            mw = Descriptors.MolWt(mol)
            if 100 <= mw <= 1500:
                neg_smiles_list.append(std_smiles)
                neg_inchikey_list.append(std_inchikey)
            else:
                print("Molecular weight out of range (100-1500):", std_smiles, mw)
        #if i > 10000:
        #    break

ik2smi_neg = {}
for smi, ik in zip(neg_smiles_list, neg_inchikey_list):
    ik2smi_neg[ik] = smi

with open(os.path.join(root, "..", "results", "00_background.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["smiles", "inchikey"])
    for ik, smi in tqdm(ik2smi_neg.items(), desc="Writing background compounds"):
        writer.writerow([smi, ik])