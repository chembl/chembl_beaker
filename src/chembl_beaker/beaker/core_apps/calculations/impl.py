from beaker.utils.io import _parseMolData
from chembl_structure_pipeline.standardizer import parse_molblock
from rdkit import Chem
import io


def _getSDFString(mols):
    sio = io.StringIO()
    _create_sdf(sio, mols)
    return sio.getvalue()


def _create_sdf(f, mols):
    for mol, props in mols:
        m = Chem.MolToMolBlock(mol, kekulize=True)
        if props:
            m += props
        f.write(m + '\n$$$$\n')

# ----------------------------------------------------------------------------------------------------------------------


def remove_hs_from_mol(m):
    indices = []
    for atom in m.GetAtoms():
        if atom.GetAtomicNum() == 1 and not atom.GetIsotope():
            bnd = atom.GetBonds()[0]
            if not (bnd.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)) and \
                    not (bnd.HasProp("_MolFileBondStereo") and bnd.GetUnsignedProp("_MolFileBondStereo") in (1, 6)):
                indices.append(atom.GetIdx())

    mol = Chem.RWMol(m)
    for index in sorted(indices, reverse=True):
        mol.RemoveAtom(index)
    return mol


def _removeHs(data):
    mols = _parseMolData(data, loadMol=False, useRDKitChemistry=False)
    ms = []
    for molblock in mols:
        props = molblock.split("M  END")[1].strip()
        props =  props if len(props) > 1 else None
        mol = parse_molblock(molblock, useRDKitChemistry=False)
        Chem.FastFindRings(mol)
        mol.UpdatePropertyCache(strict=False)
        mol = remove_hs_from_mol(mol)
        ms.append((mol, props))
    return _getSDFString(ms)
