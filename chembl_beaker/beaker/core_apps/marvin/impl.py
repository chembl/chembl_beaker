__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_beaker.beaker.core_apps.marvin.MarvinJSONEncoder import MolToMarvin, MarvinToMol

#-----------------------------------------------------------------------------------------------------------------------

def _clean2D(mrv):
    block = MarvinToMol(mrv)
    mol = Chem.MolFromMolBlock(block)
    if not mol:
        print "No mol for block:\n %s" % block
        return mrv
    AllChem.Compute2DCoords(mol, bondLength = 0.8)
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

def _molExport(structure, input_f, output_f):

    if input_f == 'mrv':
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

    elif output_f == 'sdf':
        out_structure = Chem.MolToMolBlock(mol) + '\n$$$$\n'

    elif output_f == 'mrv':
        out_structure = MolToMarvin(Chem.MolToMolBlock(mol))

    return {"structure": out_structure, "format":output_f, "contentUrl":"","contentBaseUrl":""}

#-----------------------------------------------------------------------------------------------------------------------
