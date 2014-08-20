__author__ = 'mnowotka'

from lxml.etree import _ElementTree, Element, tostring
from lxml import objectify
from lxml.objectify import ObjectifiedElement, StringElement
import json
from chembl_beaker import __version__ as version
from rdkit import Chem
from StringIO import StringIO
import os

MOL_MARVIN_SCALE = 1.8666666517333335

molBondTypes = {
"SINGLE" : 1,
"DOUBLE" : 2,
"TRIPLE" : 3,
"AROMATIC" : 4,
"UNSPECIFIED" : 8,
}

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
                            'elementType': obj.get('elementType').split()}

                else:
                    return {"atomID": ids,
                            "x": [float(x) for x in obj.get('x3').split()],
                            'y': [float(y) for y in obj.get('y3').split()],
                            'z': [float(y) for y in obj.get('z3').split()],
                            'elementType': obj.get('elementType').split()}

            if obj.tag == 'bond':
                return {"atomRefs" : obj.get('atomRefs2').split(), "order": int(obj.get('order'))}

            if obj.tag == 'atom':
                if obj.get('x2', False):
                    return {"id" : obj.get('id'), "elementType": obj.get('elementType'), "x": float(obj.get('x2')), "y": float(obj.get('y2')), "z": 0.0}
                else:
                    return {"id" : obj.get('id'), "elementType": obj.get('elementType'), "x": float(obj.get('x3')), "y": float(obj.get('y3')), "z": float(obj.get('z3'))}

        return json.JSONEncoder.default(self, obj)

#-----------------------------------------------------------------------------------------------------------------------

def _dict2Mol(obj, scale = 1.0):
    buffer = StringIO()
    buffer.write('\nConverted by chembl_beaker ver. %s\n\n' % version)
    atoms = obj['atoms']
    if type(atoms) == dict:
        num_atoms = len(atoms['atomID'])
    else:
        num_atoms = len(atoms)
    zeros = [0] * num_atoms
    buffer.write('{:3}{:3}{:3}{:3}{:3}{:3}          {} V2000\n'.format(
                                                            num_atoms, len(obj['bonds']), 0, 0, 0, 0, 1))

    if type(atoms) == dict:
        for i, atom in enumerate(atoms['elementType']):
            buffer.write('{:10,.4f}{:10,.4f}{:10,.4f} {:<2}{:3}{:3}\n'.format(
                          obj['atoms'].get('x', zeros)[i] / scale,
                          obj['atoms'].get('y', zeros)[i] / scale,
                          obj['atoms'].get('z', zeros)[i] / scale,
                          atom, 0, 0 ))
    else:
        for i, atom in enumerate(atoms):
            buffer.write('{:10,.4f}{:10,.4f}{:10,.4f} {:<2}{:3}{:3}\n'.format(
                          atom.get('x', 0) / scale,
                          atom.get('y', 0) / scale,
                          atom.get('z', 0) / scale,
                          atom.get('elementType'), 0, 0 ))

    for bond in obj['bonds']:
        first_atom = None
        second_atom = None
        first_id, second_id = bond['atomRefs']
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
        buffer.write('{:3}{:3}{:3}{:3}\n'.format(first_atom, second_atom, bond['order'], 0))
    buffer.write('M  END\n')
    return buffer.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

def _jsonToMol(obj, scale = 1.0):
    return '$$$$\n'.join([_dict2Mol(x.values()[0], scale) for x in obj])

#-----------------------------------------------------------------------------------------------------------------------

def _molToJson(mol, scale = 1.0):
    atoms = {"atomID":[], "elementType":[], "x":[], "y":[]}
    conformer = mol.GetConformer()
    is3D = conformer.Is3D()
    if is3D:
        atoms['z'] = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        point = conformer.GetAtomPosition(idx)
        atoms['x'].append(point.x * scale)
        atoms['y'].append(point.y * scale)
        if is3D:
            atoms['z'].append(point.y * scale)
        atoms["atomID"].append("a%s" % (idx + 1))
        atoms["elementType"].append(atom.GetSymbol())
    bonds = []
    for bond in mol.GetBonds():
        atomRefs = ["a%s" % (bond.GetBeginAtomIdx() + 1), "a%s" % (bond.GetEndAtomIdx() + 1)]
        order = molBondTypes[str(bond.GetBondType())]
        bonds.append({"atomRefs": atomRefs, "order":order})
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
            kwargs = {}
            kwargs['atomID'] = ' '.join(data['atomID'])
            kwargs['elementType'] = ' '.join(data['elementType'])
            if 'z' not in data:
                kwargs['x2'] = ' '.join("%.15f" % x for x in data['x'])
                kwargs['y2'] = ' '.join("%.15f" % y for y in data['y'])
            else:
                kwargs['x3'] = ' '.join("%.15f" % x for x in data['x'])
                kwargs['y3'] = ' '.join("%.15f" % y for y in data['y'])
                kwargs['z3'] = ' '.join("%.15f" % z for z in data['z'])

            element = Element('atomArray', **kwargs)
        elif name == 'bonds':
            element = Element('bondArray')
            for bond in data:
                element.append(_dictToEtree(bond, 'bond', depth + 1))
    elif depth == 5:
        kwargs = {}
        kwargs['atomRefs2'] = ' '.join(data['atomRefs'])
        kwargs['order'] = str(data['order'])
        element = Element('bond', **kwargs)
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

def SDFToMarvin(mol):
    mols = []
    if mol.endswith('.sdf') and os.path.exists(mol):
        suppl = Chem.SDMolSupplier('mol')
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