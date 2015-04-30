__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import MCS
from chembl_beaker.beaker.utils.io import _parseMolData

#-----------------------------------------------------------------------------------------------------------------------

def _mcs(data, asSmiles, atomCompare, bondCompare, threshold, ringMatchesRingOnly, completeRingsOnly, sanitize=True,
         removeHs=True, strictParsing=True, isomericSmiles=False, canonical=True, kekuleSmiles=False):
    ms = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    if not ms:
        return
    if len(ms) == 1:
        if asSmiles:
            print 'SMARTS'
            return Chem.MolToSmiles(ms[0])
        else:
            print 'SMILES'
            return Chem.MolToSmarts(ms[0])

    if threshold:
        threshold=float(threshold)
    try:
        mcs = MCS.FindMCS(ms,
                          atomCompare=atomCompare,
                          bondCompare=bondCompare,
                          ringMatchesRingOnly=ringMatchesRingOnly,
                          completeRingsOnly=completeRingsOnly,
                          threshold=threshold)
    except TypeError:
        mcs = MCS.FindMCS(ms,
                          atomCompare=atomCompare,
                          bondCompare=bondCompare,
                          ringMatchesRingOnly=ringMatchesRingOnly,
                          completeRingsOnly=completeRingsOnly)
    res = mcs.smarts
    if asSmiles:
        p = Chem.MolFromSmarts(res)
        for m in ms:
            if m.HasSubstructMatch(p):
                match = m.GetSubstructMatch(p)
                res = Chem.MolFragmentToSmiles(m, atomsToUse=match, isomericSmiles=isomericSmiles, canonical=canonical,
                kekuleSmiles=kekuleSmiles)
                break
    return res

#-----------------------------------------------------------------------------------------------------------------------