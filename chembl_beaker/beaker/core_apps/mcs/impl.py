__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
try:
    from rdkit.Chem import rdFMCS as MCS
except:
    from rdkit.Chem import MCS
from chembl_beaker.beaker.utils.io import _parseMolData


# ----------------------------------------------------------------------------------------------------------------------


def _mcs(data, asSmiles, atomCompare, bondCompare, threshold, ringMatchesRingOnly, completeRingsOnly, sanitize=True,
         removeHs=True, strictParsing=True, isomericSmiles=False, canonical=True, kekuleSmiles=False):
    ms = _parseMolData(data, sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
    if not ms:
        return
    if len(ms) == 1:
        if asSmiles:
            return Chem.MolToSmiles(ms[0])
        else:
            return Chem.MolToSmarts(ms[0])

    if threshold:
        threshold = float(threshold)
    try:
        mcs = MCS.FindMCS(ms,
                          atomCompare=atomCompare,
                          bondCompare=bondCompare,
                          ringMatchesRingOnly=ringMatchesRingOnly,
                          completeRingsOnly=completeRingsOnly,
                          threshold=threshold)
    except TypeError:
        ac = MCS.AtomCompare.CompareAny
        if hasattr(MCS.AtomCompare, atomCompare):
            ac = getattr(MCS.AtomCompare, atomCompare)
        bc = MCS.BondCompare.CompareOrder
        if hasattr(MCS.BondCompare, bondCompare):
            bc = getattr(MCS.BondCompare, bondCompare)
        th = 1.0
        if threshold:
            th = threshold
        mcs = MCS.FindMCS(ms,
                          atomCompare=ac,
                          bondCompare=bc,
                          ringMatchesRingOnly=ringMatchesRingOnly,
                          completeRingsOnly=completeRingsOnly,
                          threshold=th
                          )
    if hasattr(mcs, 'smarts'):
        res = mcs.smarts
    else:
        res = mcs.smartsString
    if asSmiles:
        p = Chem.MolFromSmarts(res)
        for m in ms:
            if m.HasSubstructMatch(p):
                match = m.GetSubstructMatch(p)
                res = Chem.MolFragmentToSmiles(m, atomsToUse=match, isomericSmiles=isomericSmiles, canonical=canonical,
                                               kekuleSmiles=kekuleSmiles)
                break
    return res

# ----------------------------------------------------------------------------------------------------------------------

