__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_beaker.beaker.core_apps.marvin.MarvinJSONEncoder import MolToMarvin, MarvinToMol
from chembl_beaker.beaker.core_apps.D3Coords.impl import _2D23D

def _hydrogenize(block, hydro):
    mol = Chem.MolFromMolBlock(block)
    res = Chem.AddHs(mol) if hydro else Chem.RemoveHs(mol)
    AllChem.Compute2DCoords(res, bondLength = 0.8)
    return MolToMarvin(Chem.MolToMolBlock(res))

#-----------------------------------------------------------------------------------------------------------------------

def _clean(mrv, dim=2):
    block = MarvinToMol(mrv)
    mol = Chem.MolFromMolBlock(block)
    if not mol:
        print "No mol for block:\n %s" % block
        return mrv
    AllChem.Compute2DCoords(mol, bondLength = 0.8)
    if dim == 3:
        mol = _2D23D(mol, True)
        mol = Chem.RemoveHs(mol)
    return MolToMarvin(Chem.MolToMolBlock(mol))

#-----------------------------------------------------------------------------------------------------------------------

def _stereoInfo(mrv):
    ret = {"headers":
               {"tetraHedral":
                    {"name":"tetraHedral","type":"COMPLEX","source":"CALCULATOR"},
                "doubleBond":
                    {"name":"doubleBond","type":"COMPLEX","source":"CALCULATOR"}
               },
           "tetraHedral":[],
           "doubleBond":[]
    }

    block = MarvinToMol(mrv)
    mol = Chem.MolFromMolBlock(block)
    if not mol:
        print "No mol for block:\n %s" % block
        return ret
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    for atom in mol.GetAtoms():
        stereo = str(atom.GetChiralTag())
        atomIndex = atom.GetIdx()
        if str(atom.GetChiralTag()) != "CHI_UNSPECIFIED":
            if stereo == "CHI_TETRAHEDRAL_CW":
                chirality = "R"
            elif stereo == "CHI_TETRAHEDRAL_CCW":
                chirality = "S"
            else:
                chirality = "R/S"
            ret["tetraHedral"].append({"atomIndex":atomIndex,"chirality":chirality})
        elif atom.HasProp('_ChiralityPossible'):
            chirality = "R/S"
            ret["tetraHedral"].append({"atomIndex":atomIndex,"chirality":chirality})


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
            ret["doubleBond"].append({"atom1Index": atom1Index,"atom2Index":atom2Index,"bondIndex":bondIndex,
                                                                                                "cistrans":cistrans})
    return ret

#-----------------------------------------------------------------------------------------------------------------------

def _molExport(structure, **kwargs):

    input_f = kwargs.get('input', None)
    output_f = kwargs.get('output', None)

    if not input_f:
        mol = _autoDetect(str(structure))
    elif input_f == 'mrv':
        mol = Chem.MolFromMolBlock(MarvinToMol(structure))
    elif input_f == 'smiles':
        mol = Chem.MolFromSmiles(str(structure))
    elif input_f == 'inchi':
        mol = Chem.MolFromInchi(str(structure))
    else:
        mol = Chem.MolFromMolBlock(structure)

    if not mol.GetNumConformers() or mol.GetConformer().Is3D():
        AllChem.Compute2DCoords(mol, bondLength = 0.8)

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

    return {"structure": out_structure, "format":output_f, "contentUrl":"","contentBaseUrl":""}

#-----------------------------------------------------------------------------------------------------------------------

def _autoDetect(structure):

    if Chem.MolFromSmiles(structure):
        return Chem.MolFromSmiles(structure.strip())

    if Chem.MolFromMolBlock(structure):
        return Chem.MolFromMolBlock(structure)

    if Chem.inchi.MolFromInchi(structure.strip(), True, True):
        return Chem.inchi.MolFromInchi(structure.strip(), True, True)

    if Chem.MolFromSmarts(structure):
        return Chem.MolFromSmarts(structure)

    if Chem.MolFromMol2Block(structure):
        return Chem.MolFromMol2Block(structure)

    if Chem.MolFromPDBBlock(structure):
        return Chem.MolFromPDBBlock(structure)

    if Chem.MolFromTPLBlock(structure):
        return Chem.MolFromTPLBlock(structure)

    return None

#-----------------------------------------------------------------------------------------------------------------------
