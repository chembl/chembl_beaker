__author__ = 'mnowotka'

from bottle import Bottle, route, run, get, post, request, response, default_app
from chembl_beaker import settings
from chembl_beaker import __version__ as version
import base64
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import StringIO
import os
from subprocess import PIPE, Popen
import tempfile
import json

app = Bottle()

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/status')
def status():
    return "This is ChEMBL beaker, version %s" % version

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2smiles/<ctab>')
def ctab2smiles(ctab):
    data = base64.b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    sio = StringIO.StringIO()
    w = Chem.SmilesWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2smiles')
def ctab2smiles():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    mols = [x for x in suppl]
    sio = StringIO.StringIO()
    w = Chem.SmilesWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/smiles2ctab/<smiles>')
def smiles2ctab(smiles):
    data = base64.b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    suppl = Chem.SmilesMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    for m in mols:
        AllChem.Compute2DCoords(m)
       
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/smiles2ctab')
def smiles2ctab():
    data = request.body.getvalue()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    suppl = Chem.SmilesMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    for m in mols:
        AllChem.Compute2DCoords(m)
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/inchi2ctab/<inchi>')
def inchi2ctab(inchi):
    inchis = base64.b64decode(inchi)
    mols = [Chem.MolFromInchi(inch) for inch in inchis.split()]
    for m in mols:
        AllChem.Compute2DCoords(m)
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/inchi2ctab')
def inchi2ctab():
    inchis = request.body.getvalue()
