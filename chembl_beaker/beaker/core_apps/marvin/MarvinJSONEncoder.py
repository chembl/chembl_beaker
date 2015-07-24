__author__ = 'mnowotka'

from lxml.etree import _ElementTree, Element, tostring
from lxml import objectify
from lxml.objectify import ObjectifiedElement, StringElement
import json
from chembl_beaker import __version__ as version
from datetime import datetime, date
import rdkit
from rdkit import Chem
from rdkit.Chem.rdchem import GetPeriodicTable
from StringIO import StringIO
import os
import getpass

USERNAME = ''
try:
    USERNAME = getpass.getuser()
except:
    pass

MOL_MARVIN_SCALE = 1.8666666517333335

molBondTypes = {
"SINGLE" : 1,
"DOUBLE" : 2,
"TRIPLE" : 3,
"AROMATIC" : 4,
"UNSPECIFIED" : 8,
}

charges = {
    0:0,
    1:3,
    2:2,
    3:1,
    4:0,
    5:-1,
    6:-2,
    7:-3,
}

PERIODIC_TABLE = GetPeriodicTable()

#-----------------------------------------------------------------------------------------------------------------------

class MarvinJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if type(obj) == _ElementTree:
            return obj.getroot()
        if type(obj) == ObjectifiedElement:
            if obj.tag == 'cml':
                return  obj.find('MDocument')
            if obj.tag == 'MDocument':
                return [ el for el in obj.iterchildren(tag='MChemicalStruct') ]
            if obj.tag == 'MChemicalStruct':
                return obj.find('molecule')
            if obj.tag == 'molecule':
                return {obj.get("molID"): {"atoms" : obj.find("atomArray"), "bonds": obj.find("bondArray")}}
            if obj.tag == 'bondArray':
                return [ bond for bond in obj.iterchildren(tag='bond') ]
            if obj.tag == 'atomArray':
                return [ atom for atom in obj.iterchildren(tag='atom') ]
            if obj.tag == 'bond':
                return {"atomRefs" : obj.get('atomRefs2').split(),
                        "order": int(obj.get('order')),
                        "stereo": obj.find("bondStereo")}
            if obj.tag == 'bondStereo':
                pass
        if type(obj) == StringElement:
            if obj.tag == 'atomArray':
                ids = obj.get('atomID').split()
                num_atoms = len(ids)
                zeros = [0] * num_atoms

                if obj.get('x2', False):
                    return {"atomID": ids,
                            "x": [float(x) for x in obj.get('x2').split()],
                            'y': [float(y) for y in obj.get('y2').split()],
                            'z': zeros,
                            'dim': 2,
                            'elementType': obj.get('elementType').split()}

                else:
                    return {"atomID": ids,
                            "x": [float(x) for x in obj.get('x3').split()],
                            'y': [float(y) for y in obj.get('y3').split()],
                            'z': [float(y) for y in obj.get('z3').split()],
                            'dim':3,
                            'elementType': obj.get('elementType').split()}

            if obj.tag == 'bond':
                return {"atomRefs" : obj.get('atomRefs2').split(), "order": int(obj.get('order'))}

            if obj.tag == 'atom':
                if obj.get('x2', False):
                    return {"id" : obj.get('id'),
                            "elementType": obj.get('elementType'),
                            "x": float(obj.get('x2')),
                            "y": float(obj.get('y2')),
                            "z": 0.0,
                            'dim':2,
                            'formalCharge': int(obj.get('formalCharge', 0)),
                            'isotope':int(obj.get('isotope', 0))
                    }
                else:
                    return {"id" : obj.get('id'),
                            "elementType": obj.get('elementType'),
                            "x": float(obj.get('x3')),
                            "y": float(obj.get('y3')),
                            "z": float(obj.get('z3')),
                            'dim':3,
                            'formalCharge': int(obj.get('formalCharge', 0)),
                            'isotope':int(obj.get('isotope', 0))}

            if obj.tag == 'bondStereo':
                convention = obj.get('convention')
                if convention and convention.upper() == 'MDL':
                    return int(obj.get('conventionValue', 0))
                if not obj.text:
                    return 0
                t = obj.text.strip().upper()
                if t == 'W':
                    return 1
                elif t == 'H':
                    return 6
                else:
                    return 0

            if obj.tag == 'bondArray':
                return []

        return json.JSONEncoder.default(self, obj)

#-----------------------------------------------------------------------------------------------------------------------

def _dict2Mol(obj, scale = 1.0):
    buffer = StringIO()
    data = date.today().strftime("%m%d%y%H%M")

    atoms = obj['atoms']
    if type(atoms) == dict:
        num_atoms = len(atoms['atomID'])
    else:
        num_atoms = len(atoms)
    if type(atoms) == dict:
        dim = str(atoms.get('dim', 2)) + 'D'
    else:
        dim = str(atoms[0].get('dim',2)) + 'D'

    buffer.write('\n')
    buffer.write('{:2.2}{:^8.8}{:10.10}{:2.2}{:2.2}{:10.10}{:12.12}{:6.6}\n'.format(
                                                                    USERNAME, 'beaker', data, dim, '', '', '', ''))
    buffer.write('\n')

    zeros = [0] * num_atoms
    buffer.write('{:3}{:3}{:3}{:3}{:3}{:3}          {} V2000\n'.format(
                                                            num_atoms, len(obj['bonds']), 0, 0, 0, 0, 1))

    if type(atoms) == dict:
        for i, atom in enumerate(atoms['elementType']):
            buffer.write('{:10,.4f}{:10,.4f}{:10,.4f} {:<2}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}\n'.format(
                          atoms.get('x', zeros)[i] / scale,
                          atoms.get('y', zeros)[i] / scale,
                          atoms.get('z', zeros)[i] / scale,
                          atom,
                          0 if not atoms.get('isotope') else atoms.get('isotope')[i] -
                                                             int(PERIODIC_TABLE.GetAtomicWeight(str(atom))),
                          charges.keys()[charges.values().index(atoms.get('formalCharge', zeros)[i])],
                                                                                                0,0,0,0,0,0,0,0,0,0))
    else:
        for i, atom in enumerate(atoms):
            elementType = atom.get('elementType')
            buffer.write('{:10,.4f}{:10,.4f}{:10,.4f} {:<2}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}{:3}\n'.format(
                          atom.get('x', 0) / scale,
                          atom.get('y', 0) / scale,
                          atom.get('z', 0) / scale,
                          elementType,
                          0 if not atom.get('isotope') else atom.get('isotope') -
                                                            int(PERIODIC_TABLE.GetAtomicWeight(str(elementType))),
                          charges.keys()[charges.values().index(atom.get('formalCharge', 0))],
                                                                                                0,0,0,0,0,0,0,0,0,0 ))

    for bond in obj['bonds']:
        first_atom = None
        second_atom = None
        first_id, second_id = bond['atomRefs']
        stereo = bond.get('stereo', 0)
        if type(atoms) == dict:
            for i, atomID in enumerate(atoms['atomID']):
                if first_id == atomID:
                    first_atom = i + 1
                if second_id == atomID:
                    second_atom = i + 1
                if first_atom is not None and second_atom is not None:
                    break
        else:
            for i, atom in enumerate(atoms):
                atomID = atom.get('id')
                if first_id == atomID:
                    first_atom = i + 1
                if second_id == atomID:
                    second_atom = i + 1
                if first_atom is not None and second_atom is not None:
                    break
        buffer.write('{:3}{:3}{:3}{:3}{:3}{:3}{:3}\n'.format(first_atom, second_atom, bond['order'], stereo, 0, 0, 0))
    buffer.write('M  END\n')
    return buffer.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _jsonToMol(obj, scale = 1.0):
    return '$$$$\n'.join([_dict2Mol(x.values()[0], scale) for x in obj])

#-----------------------------------------------------------------------------------------------------------------------

def _molToJson(mol, scale = 1.0):
    atoms = []
    conformer = mol.GetConformer()
    is3D = conformer.Is3D()
    for atom in mol.GetAtoms():
        atom_data = {}
        idx = atom.GetIdx()
        point = conformer.GetAtomPosition(idx)
        x = point.x * scale
        y = point.y * scale
        if is3D:
            z = point.y * scale
            atom_data['x3'] = str(x)
            atom_data['y3'] = str(y)
            atom_data['z3'] = str(z)
        else:
            atom_data['x2'] = str(x)
            atom_data['y2'] = str(y)
        atom_data['id'] = "a%s" % (idx + 1)
        atom_data["elementType"] = atom.GetSymbol()
        formalCharge = atom.GetFormalCharge()
        if formalCharge:
            atom_data["formalCharge"] = str(formalCharge)
        isotope = atom.GetIsotope()
        if isotope:
            atom_data["isotope"] = str(isotope)
        atoms.append(atom_data)
    bonds = []
    for bond in mol.GetBonds():
        atomRefs = ["a%s" % (bond.GetBeginAtomIdx() + 1), "a%s" % (bond.GetEndAtomIdx() + 1)]
        order = molBondTypes[str(bond.GetBondType())]
        stereo = bond.GetBondDir()
        bonds.append({"atomRefs": atomRefs, "order":order, "stereo":stereo})
    return {"atoms":atoms, "bonds":bonds}

#-----------------------------------------------------------------------------------------------------------------------

def _molsToJson(mols, scale = 1.0):
    ret = []
    for idx, mol in enumerate(mols):
        ret.append({"m%s" % (idx + 1): _molToJson(mol, scale)})
    return ret

#-----------------------------------------------------------------------------------------------------------------------

def _dictToEtree(data, name=None, depth=0):
    element = None
    if depth == 0:
        element = Element('cml')
        element.append(_dictToEtree(data, name, depth + 1))
    elif depth == 1:
        element = Element('MDocument')
        element.append(_dictToEtree(data, name, depth + 1))
    elif depth == 2:
        element = Element('MChemicalStruct')
        for mol in data:
            element.append(_dictToEtree(mol, name, depth + 1))
    elif depth == 3:
        molID = data.keys()[0]
        val = data.values()[0]
        element = Element('molecule', molID=molID)
        element.append(_dictToEtree(val['atoms'], 'atoms', depth + 1))
        element.append(_dictToEtree(val['bonds'], 'bonds', depth + 1))
    elif depth == 4:
        if name == 'atoms':
            element = Element('atomArray')
            for atom in data:
                element.append(_dictToEtree(atom, 'atom', depth + 1))
        elif name == 'bonds':
            element = Element('bondArray')
            for bond in data:
                element.append(_dictToEtree(bond, 'bond', depth + 1))
    elif depth == 5:
        if name == 'bond':
            kwargs = {}
            kwargs['atomRefs2'] = ' '.join(data['atomRefs'])
            kwargs['order'] = str(data['order'])
            element = Element('bond', **kwargs)
            stereo = data.get('stereo')
            if stereo:
                element.append(_dictToEtree(stereo, 'bondStereo', depth + 1))
        elif name == 'atom':
            element = Element('atom', **data)
    elif depth == 6:
        element = Element('bondStereo')
        if data == rdkit.Chem.rdchem.BondDir.BEGINWEDGE:
            element.text = 'W'
        elif data == rdkit.Chem.rdchem.BondDir.BEGINDASH:
            element.text = 'H'
        elif data == rdkit.Chem.rdchem.BondDir.UNKNOWN:
            element.attrib["convention"] = "MDL"
            element.attrib["conventionValue"] = "4"
    return element

#-----------------------------------------------------------------------------------------------------------------------

def _dataToXml(data, options=None):
    options = options or {}
    return tostring(_dictToEtree(data, options), xml_declaration=False, encoding='utf-8')

#-----------------------------------------------------------------------------------------------------------------------

def MolToMarvin(mol):
    if mol.endswith('.mol') and os.path.exists(mol):
        mol = Chem.MolFromMolFile(mol, False, False, False)
    else:
        mol = Chem.MolFromMolBlock(mol, False, False, False)
    js = _molsToJson([mol], MOL_MARVIN_SCALE)
    return _dataToXml(js)

#-----------------------------------------------------------------------------------------------------------------------

def SDFToMarvin(mol, removeHs=True):
    mols = []
    if mol.endswith('.sdf') and os.path.exists(mol):
        suppl = Chem.SDMolSupplier('mol', removeHs=removeHs)
        mols = [mol for mol in suppl]

    js = _molsToJson(mols, MOL_MARVIN_SCALE)
    return _dataToXml(js)

#-----------------------------------------------------------------------------------------------------------------------

def MarvinToMol(marvin):
    if marvin.endswith('.mrv') and os.path.exists(marvin):
        f = open(marvin)
    else:
        f = StringIO(marvin)
    tree = objectify.parse(f)
    js = MarvinJSONEncoder().encode(tree)
    return _jsonToMol(json.loads(js), MOL_MARVIN_SCALE)

#-----------------------------------------------------------------------------------------------------------------------