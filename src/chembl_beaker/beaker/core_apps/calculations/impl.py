from beaker.utils.io import _parseMolData
from beaker.core_apps.D2Coords.impl import _check3Dcoords
from chembl_structure_pipeline.standardizer import remove_hs_from_mol
from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit import Chem
import io


def _getSDFString(mols):
    sio = io.StringIO()
    _create_sdf(sio, mols)
    return sio.getvalue()


def _create_sdf(f, mols):
    for mol, props in mols:
        m = Chem.MolToMolBlock(mol)
        if props:
            m += props
        f.write(m + '\n$$$$\n')

# ----------------------------------------------------------------------------------------------------------------------

def _removeHs(data):
    mols = _parseMolData(data, loadMol=False, useRDKitChemistry=False)
    ms = []
    for molblock in mols:
        mol = parse_molblock(molblock, useRDKitChemistry=False)
        props = molblock.split("M  END")[1].strip()
        props =  props if len(props) > 1 else None
        Chem.FastFindRings(mol)
        if mol.NeedsUpdatePropertyCache():
            mol.UpdatePropertyCache(strict=False)
        if _check3Dcoords(mol):
            Chem.AssignStereochemistryFrom3D(mol)
        else:
            Chem.AssignAtomChiralTagsFromStructure(mol)
        ms.append((remove_hs_from_mol(mol), props))
    return _getSDFString(ms)
