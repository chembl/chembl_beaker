__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_beaker.beaker.core_apps.marvin.MarvinJSONEncoder import MolToMarvin, MarvinToMol
from chembl_beaker.beaker.core_apps.D3Coords.impl import _2D23D
from chembl_beaker.beaker.utils.functional import _apply, _call
import chembl_beaker.beaker.utils.chemical_transformation as ct


def _hydrogenize(block, hydro):
    mol = Chem.MolFromMolBlock(block, sanitize=False)
    _call([mol], 'UpdatePropertyCache', strict=False)
    _apply([mol], ct._sssr)
    res = Chem.AddHs(mol, addCoords=True) if hydro else Chem.RemoveHs(mol, sanitize=False)
    return MolToMarvin(Chem.MolToMolBlock(res))


# ----------------------------------------------------------------------------------------------------------------------


def _clean(mrv, dim=2):
    block = MarvinToMol(mrv)
    mol = Chem.MolFromMolBlock(block, sanitize=False)
    _call([mol], 'UpdatePropertyCache', strict=False)
    _apply([mol], ct._sssr)    
    if not mol:
        print "No mol for block:\n %s" % block
        return mrv
    AllChem.Compute2DCoords(mol, bondLength=0.8)
    if dim == 3:
        mol = _2D23D(mol, True)
        mol = Chem.RemoveHs(mol)
    return MolToMarvin(Chem.MolToMolBlock(mol))


# ----------------------------------------------------------------------------------------------------------------------


def _stereoInfo(mrv):
    ret = {"headers":
               {"tetraHedral":
                    {"name": "tetraHedral", "type": "COMPLEX", "source": "CALCULATOR"},
                "doubleBond":
                    {"name": "doubleBond", "type": "COMPLEX", "source": "CALCULATOR"}
                },
           "tetraHedral": [],
           "doubleBond": []
           }

    block = MarvinToMol(mrv)
    mol = Chem.MolFromMolBlock(block)
    if not mol:
        print "No mol for block:\n %s" % block
        return ret
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    for atom in mol.GetAtoms():
        atomIndex = atom.GetIdx()
        if atom.HasProp('_CIPCode'):
            ret["tetraHedral"].append({"atomIndex": atomIndex, "chirality": atom.GetProp('_CIPCode')})

    for bond in mol.GetBonds():
        stereo = str(bond.GetStereo())
        if stereo != "STEREONONE":
            atom1Index = bond.GetBeginAtomIdx()
            atom2Index = bond.GetEndAtomIdx()
            bondIndex = bond.GetIdx()
            if stereo == "STEREOANY":
                cistrans = "E/Z"
            elif stereo == "STEREOZ":
                cistrans = "Z"
            elif stereo == "STEREOE":
                cistrans = "E"
            ret["doubleBond"].append({"atom1Index": atom1Index, "atom2Index": atom2Index, "bondIndex": bondIndex,
                                      "cistrans": cistrans})
    return ret


# ----------------------------------------------------------------------------------------------------------------------

def _molExport(structure, **kwargs):
    input_f = kwargs.get('input', None)
    output_f = kwargs.get('output', None)

    if not input_f:
        mol = _autoDetect(str(structure))
    elif input_f == 'mrv':
        mol = Chem.MolFromMolBlock(MarvinToMol(structure), sanitize=False)
    elif input_f == 'smiles':
        mol = Chem.MolFromSmiles(str(structure), sanitize=False)
    elif input_f == 'inchi':
        mol = Chem.MolFromInchi(str(structure))
    else:
        mol = Chem.MolFromMolBlock(structure, sanitize=False)

    _call([mol], 'UpdatePropertyCache', strict=False)
    _apply([mol], ct._sssr)

    if not mol.GetNumConformers() or mol.GetConformer().Is3D():
        AllChem.Compute2DCoords(mol, bondLength=0.8)

    if output_f == 'smiles':
        out_structure = Chem.MolToSmiles(mol)

    elif output_f == 'inchi':
        out_structure = Chem.MolToInchi(mol)

    elif output_f == 'inchikey':
        out_structure = Chem.InchiToInchiKey(Chem.MolToInchi(mol))

    elif output_f == 'mol':
        out_structure = Chem.MolToMolBlock(mol)

    elif output_f == 'sdf':
        out_structure = Chem.MolToMolBlock(mol) + '\n$$$$\n'

    elif output_f == 'mrv':
        out_structure = MolToMarvin(Chem.MolToMolBlock(mol))

    return {"structure": out_structure, "format": output_f, "contentUrl": "", "contentBaseUrl": ""}


# ----------------------------------------------------------------------------------------------------------------------

def _autoDetect(structure):
    mol = None

    if Chem.MolFromSmiles(structure.strip(), sanitize=False):
        mol = Chem.MolFromSmiles(structure.strip(), sanitize=False)

    elif Chem.MolFromMolBlock(structure, sanitize=False):
        return Chem.MolFromMolBlock(structure, sanitize=False)

    elif Chem.INCHI_AVAILABLE and Chem.inchi.MolFromInchi(structure.strip(), sanitize=True, removeHs=True):
        mol = Chem.inchi.MolFromInchi(structure.strip(), sanitize=True, removeHs=True)

    elif hasattr(Chem, 'MolFromSmarts') and Chem.MolFromSmarts(structure):
        mol = Chem.MolFromSmarts(structure)

    elif hasattr(Chem, 'MolFromMol2Block') and Chem.MolFromMol2Block(structure):
        mol = Chem.MolFromMol2Block(structure)

    elif hasattr(Chem, 'MolFromPDBBlock') and Chem.MolFromPDBBlock(structure):
        mol = Chem.MolFromPDBBlock(structure)

    elif hasattr(Chem, 'MolFromTPLBlock') and Chem.MolFromTPLBlock(structure):
        mol = Chem.MolFromTPLBlock(structure)

    else:
        try:
            mol = Chem.MolFromMolBlock(MarvinToMol(structure))
        except:
            pass

    _call([mol], 'UpdatePropertyCache', strict=False)
    _apply([mol], ct._sssr)

    return mol

# ----------------------------------------------------------------------------------------------------------------------
