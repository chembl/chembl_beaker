from beaker.utils.functional import _apply
from beaker.utils.io import _parseMolData, _getSDFString
from beaker.core_apps.D2Coords.impl import _check3Dcoords
from chembl_structure_pipeline.standardizer import remove_hs_from_mol
from rdkit import Chem

# ----------------------------------------------------------------------------------------------------------------------


def _removeHs(data):
    mols = _parseMolData(data, loadMol=True, useRDKitChemistry=False)
    for mol in mols:
        Chem.FastFindRings(mol)
        if mol.NeedsUpdatePropertyCache():
            mol.UpdatePropertyCache(strict=False)
        if _check3Dcoords(mol):
            Chem.AssignStereochemistryFrom3D(mol)
        else:
            Chem.AssignAtomChiralTagsFromStructure(mol)
    mols = _apply(mols, remove_hs_from_mol)
    return _getSDFString(mols)
