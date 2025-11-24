import os
import csv

from FPSim2.io import create_db_file

root = os.path.join(os.path.dirname(__file__))

mols = []
with open(os.path.join(root, '..', 'results', '00_abx_compounds.csv'), 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for i, r in enumerate(reader):
        smi = r[0]
        mols.append([smi, i+1])

checkpoints_dir = os.path.join(root, '..', '..', '..', 'checkpoints')

def morgan_fp():
    create_db_file(
        mols_source=mols,
        filename=os.path.join(checkpoints_dir, 'morgan_2048.h5'),
        mol_format='smiles',
        fp_type='Morgan',
        fp_params={'radius': 2, 'fpSize': 2048}
    )

def maccs_fp():
    create_db_file(
        mols_source=mols,
        filename=os.path.join(checkpoints_dir, 'maccs_166.h5'),
        mol_format='smiles',
        fp_type='MACCSKeys'
    )

def rdkit_fp():
    create_db_file(
        mols_source=mols,
        filename=os.path.join(checkpoints_dir, 'rdkit_2048.h5'),
        mol_format='smiles',
        fp_type='RDKit',
        fp_params={'fpSize': 2048}
    )


def atompair_fp():
    create_db_file(
        mols_source=mols,
        filename=os.path.join(checkpoints_dir, 'atompair_2048.h5'),
        mol_format='smiles',
        fp_type='AtomPair',
        fp_params={'fpSize': 2048}
    )

def pattern_fp():
    create_db_file(
        mols_source=mols,
        filename=os.path.join(checkpoints_dir, 'pattern_2048.h5'),
        mol_format='smiles',
        fp_type='RDKitPattern',
        fp_params={'fpSize': 2048}
    )

morgan_fp()
maccs_fp()
rdkit_fp()
atompair_fp()
pattern_fp()
