__author__ = 'efelix'

from rdkit import Chem
from beaker.utils.functional import _apply, _call
from beaker.utils.io import _parseMolData
import beaker.utils.chemical_transformation as ct
import copy
import json
import os.path as path

data_dir =  path.abspath(path.join(__file__ , "../../../../data"))

with open('{}/chembl_26_alerts.json'.format(data_dir), 'r') as alerts_file:
    alerts = json.load(alerts_file)

for alert in alerts:
    alert['rdmol'] = Chem.MolFromSmarts(alert['smarts'])

#-----------------------------------------------------------------------------------------------------------------------

def get_matches(mol):
    _call([mol], 'UpdatePropertyCache', strict=False)
    _apply([mol], ct._sssr)

    matches = []
    for alert in alerts:
        try:
            match = mol.HasSubstructMatch(alert['rdmol'])
            if match:
                malert = copy.deepcopy(alert)
                del malert['rdmol']
                matches.append(malert)
        except Exception as e:
            print(e)
    return matches

#-----------------------------------------------------------------------------------------------------------------------

def _get_alerts(data, loadMol=True, useRDKitChemistry=False):
    mols = _parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry)
    res = _apply(mols, get_matches)
    return json.dumps(res)

#-----------------------------------------------------------------------------------------------------------------------
