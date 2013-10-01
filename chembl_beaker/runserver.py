__author__ = 'mnowotka'

from bottle import route, run, get, post, request, response
from chembl_beaker import settings
import base64
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import StringIO
import os
from subprocess import PIPE, Popen
import tempfile

#-----------------------------------------------------------------------------------------------------------------------

@route('/status')
def hello():
    return "OK!"

#-----------------------------------------------------------------------------------------------------------------------

@get('/ctab2smiles/<ctab>')
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

@post('/ctab2smiles')
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

@get('/smiles2ctab/<smiles>')
def smiles2ctab(smiles):
    data = base64.b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    suppl = Chem.SmilesMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@post('/smiles2ctab')
def smiles2ctab():
    data = request.body.getvalue()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    suppl = Chem.SmilesMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@get('/inchi2ctab/<inchi>')
def inchi2ctab(inchi):
    inchis = base64.b64decode(inchi)
    mols = [Chem.MolFromInchi(inch) for inch in inchis.split()]
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@post('/inchi2ctab')
def inchi2ctab():
    inchis = request.body.getvalue()
    mols = [Chem.MolFromInchi(inch) for inch in inchis.split()]
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    for m in mols:
        w.write(m)
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@get('/ctab2inchi/<ctab>')
def ctab2inchi(ctab):
    data = base64.b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    inchi = []
    for m in mols:
        inchi.append(Chem.MolToInchi(m))
    return '\n'.join(inchi)

#-----------------------------------------------------------------------------------------------------------------------

@post('/ctab2inchi')
def ctab2inchi():
    suppl = Chem.SDMolSupplier()
    suppl.SetData(request.body.getvalue())
    mols = [x for x in suppl]
    inchi = []
    for m in mols:
        inchi.append(Chem.MolToInchi(m))
    return '\n'.join(inchi)

#-----------------------------------------------------------------------------------------------------------------------

@get('/inchi2inchiKey/<inchi>')
def inchi2inchiKey(inchi):
    inchis = base64.b64decode(inchi)
    keys = [Chem.InchiToInchiKey(inch) for inch in inchis.split()]
    return '\n'.join(keys)

#-----------------------------------------------------------------------------------------------------------------------

@post('/inchi2inchiKey')
def inchi2inchiKey():
    inchis = request.body.getvalue()
    keys = [Chem.InchiToInchiKey(inch) for inch in inchis.split()]
    return '\n'.join(keys)

#-----------------------------------------------------------------------------------------------------------------------

@get('/ctab2image/<ctab>')
@get('/ctab2image/<ctab>/<size>')
@get('/ctab2image/<ctab>/<size>/<legend>')
def ctab2image(ctab, size=200, legend=''):
    size = int(size)
    data = base64.b64decode(ctab)
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    for m in mols:
        AllChem.Compute2DCoords(m)
    image = Draw.MolsToGridImage(mols,molsPerRow=min(len(mols),4),subImgSize=(size,size),legends=[x.GetProp("_Name") or legend for x in mols])
    imageData = StringIO.StringIO()
    image.save(imageData, "PNG")
    response.content_type = 'image/png'
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@post('/ctab2image')
def ctab2image():
    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.getvalue()
    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    mols = [x for x in suppl]
    for m in mols:
        AllChem.Compute2DCoords(m)
    image = Draw.MolsToGridImage(mols,molsPerRow=min(len(mols),4),subImgSize=(size,size),legends=[x.GetProp("_Name") for x in mols])
    imageData = StringIO.StringIO()
    image.save(imageData, "PNG")
    response.content_type = 'image/png'
    return imageData.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@get('/image2ctab/<image>')
def image2ctab(image):
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    img = base64.b64decode(image)
    path = settings.OSRA_BINARIES_LOCATION['2.0.0']
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    p = Popen([path, fpath], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a, err = p.communicate(input=img)
    os.remove(fpath)
    for smiles in filter(bool,a.split('\n')):
        w.write(Chem.MolFromSmiles(smiles))
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@post('/image2ctab')
def image2ctab():
    sio = StringIO.StringIO()
    w = Chem.SDWriter(sio)
    img = request.body.read()
    path = settings.OSRA_BINARIES_LOCATION['2.0.0']
    fd, fpath = tempfile.mkstemp()
    os.write(fd, img)
    os.close(fd)
    p = Popen([path, fpath], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    a, err = p.communicate(input=img)
    os.remove(fpath)
    for smiles in filter(bool,a.split('\n')):
        w.write(Chem.MolFromSmiles(smiles))
    w.flush()
    return sio.getvalue()

#-----------------------------------------------------------------------------------------------------------------------

@get('/kekulize/<ctab>')
def kekulize(ctab):
    data = base64.b64decode(ctab)
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

@post('/kekulize')
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

@route('/sanitize/<ctab>')
def sanitize(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/atomIsInRing/<ctab>/<index>/<size>')
def sanitize(ctab, index, size):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/symmSSSR/<ctab>')
def atomIsInRing(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/SSSR/<ctab>')
def SSSR(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/addHs/<ctab>')
def addHs(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/removeHs/<ctab>')
def removeHs(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/getNumAtoms/<ctab>')
def getNumAtoms(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/logP/<ctab>')
def logP(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/TPSA/<ctab>')
def TPSA(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

@route('/molWt/<ctab>')
def molWt(ctab):
    pass

#-----------------------------------------------------------------------------------------------------------------------

run(host=settings.BOTTLE_HOST, port=settings.BOTTLE_PORT, debug=settings.DEBUG, server=settings.SERVER_MIDDLEWARE)

#-----------------------------------------------------------------------------------------------------------------------