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
    data = base64.urlsafe_b64decode(ctab)
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
    data = base64.urlsafe_b64decode(smiles)
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
    inchis = base64.urlsafe_b64decode(inchi)
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

@app.get('/ctab2inchi/<ctab>')
def ctab2inchi(ctab):
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    inchi = []
    for m in mols:
        inchi.append(Chem.MolToInchi(m))
    return '\n'.join(inchi)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2inchi')
def ctab2inchi():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    mols = [x for x in suppl]
    inchi = []
    for m in mols:
        inchi.append(Chem.MolToInchi(m))
    return '\n'.join(inchi)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/inchi2inchiKey/<inchi>')
def inchi2inchiKey(inchi):
    inchis = base64.urlsafe_b64decode(inchi)
    keys = [Chem.InchiToInchiKey(inch) for inch in inchis.split()]
    return '\n'.join(keys)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/inchi2inchiKey')
def inchi2inchiKey():
    inchis = request.body.getvalue()
    keys = [Chem.InchiToInchiKey(inch) for inch in inchis.split()]
    return '\n'.join(keys)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2image/<ctab>')
@app.get('/ctab2image/<ctab>/<size>')
@app.get('/ctab2image/<ctab>/<size>/<legend>')
def ctab2image(ctab, size=200, legend=''):
    size = int(size)
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    for m in mols:
        if m.GetConformer().Is3D():
            AllChem.Compute2DCoords(m)
    image = Draw.MolsToGridImage(mols,molsPerRow=min(len(mols),4),subImgSize=(size,size),legends=[x.GetProp("_Name") or legend for x in mols])
    imageData = StringIO.StringIO()
    image.save(imageData, "PNG")
    response.content_type = 'image/png'
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2image')
def ctab2image():
    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.getvalue()
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    for m in mols:
        if m.GetConformer().Is3D():
            AllChem.Compute2DCoords(m)
    image = Draw.MolsToGridImage(mols,molsPerRow=min(len(mols),4),subImgSize=(size,size),legends=[x.GetProp("_Name") for x in mols])
    imageData = StringIO.StringIO()
    image.save(imageData, "PNG")
    response.content_type = 'image/png'
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/image2ctab/<image>')
def image2ctab(image):
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    img = base64.urlsafe_b64decode(image)
    osras = settings.OSRA_BINARIES_LOCATION
    latest_osra = osras[max(osras.keys())]
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    p = Popen([latest_osra, fpath], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a, err = p.communicate(input=img)
    os.remove(fpath)
    for smiles in filter(bool,a.split('\n')):
        w.write(Chem.MolFromSmiles(smiles))
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/image2ctab')
def image2ctab():
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    img = request.body.read()
    osras = settings.OSRA_BINARIES_LOCATION
    latest_osra = osras[max(osras.keys())]
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    p = Popen([latest_osra, fpath], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a, err = p.communicate(input=img)
    os.remove(fpath)
    for smiles in filter(bool,a.split('\n')):
        w.write(Chem.MolFromSmiles(smiles))
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/kekulize/<ctab>')
def kekulize(ctab):
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    mols = [Chem.Kekulize(x) for x in suppl]
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/kekulize')
def kekulize():
    data = request.body.getvalue()
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    mols = [Chem.Kekulize(x) for x in suppl]
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/sanitize/<ctab>')
def sanitize(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/atomIsInRing/<ctab>/<index>/<size>')
def atomIsInRing(ctab, index, size):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/symmSSSR/<ctab>')
def symmSSSR(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/SSSR/<ctab>')
def SSSR(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/addHs/<ctab>')
def addHs(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/removeHs/<ctab>')
def removeHs(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/getNumAtoms/<ctab>')
def getNumAtoms(ctab):
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    ret = [x.GetNumAtoms() for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/getNumAtoms')
def getNumAtoms():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    ret = [x.GetNumAtoms() for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/logP/<ctab>')
def logP(ctab):
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    ret = [Descriptors.MolLogP(x) for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/logP')
def logP():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    ret = [Descriptors.MolLogP(x) for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/TPSA/<ctab>')
def TPSA(ctab):
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    ret = [Descriptors.TPSA(x) for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/TPSA')
def TPSA():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    ret = [Descriptors.TPSA(x) for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/molWt/<ctab>')
def molWt(ctab):
    data = base64.urlsafe_b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    ret = [Descriptors.MolWt(x) for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/molWt')
def molWt():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    ret = [Descriptors.MolWt(x) for x in suppl]
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    run(app=app, host=settings.BOTTLE_HOST, port=settings.BOTTLE_PORT, debug=settings.DEBUG, server=settings.SERVER_MIDDLEWARE)
else:
    application = app

#-----------------------------------------------------------------------------------------------------------------------
