__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import StringIO
from rdkit import rdBase
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolDescriptors

from chembl_beaker.beaker.utils.io import _parseMolData

#-----------------------------------------------------------------------------------------------------------------------

def _getFPSStream(f, mols, type='morgan', radius=2, n_bits=2048):
    f.write("#FPS1\n#num_bits=%s\n#software=RDKit/%s\n" % (n_bits, rdBase.rdkitVersion))
    for i, mol in enumerate(mols):
        if mol:
            idx = i
            if mol.HasProp('chembl_id'):
                idx = mol.GetProp('chembl_id')
            elif Chem.INCHI_AVAILABLE:
                try:
                    Chem.SanitizeMol(mol)
                    idx = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
                except:
                    pass
            if type == 'morgan':
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol,radius,nBits=n_bits)
            elif type == 'pair':
                fp = Pairs.GetAtomPairFingerprintAsBitVect(mol)
            elif type == 'maccs':
                fp = MACCSkeys.GenMACCSKeys(mol)
            f.write("%s\t%s\n" % (DataStructs.BitVectToFPSText(fp), idx))

#-----------------------------------------------------------------------------------------------------------------------

def _getFPSString(mols, type='morgan', radius=2, n_bits=2048):
    sio = StringIO.StringIO()
    _getFPSStream(sio, mols, type, radius, n_bits)
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _sdf2fps(sdf, type='morgan', radius=2, n_bits=2048, sanitize=True, removeHs=True, strictParsing=True):
    return _getFPSString(_parseMolData(sdf, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing),
        type, radius, n_bits)

#-----------------------------------------------------------------------------------------------------------------------
