#!/usr/bin/env python

__author__ = 'mnowotka'

import bottle
from bottle import Bottle, run, request, response
from chembl_beaker import utils
from chembl_beaker import __version__ as version
import base64
import json
import sys
from optparse import OptionParser

HTTP_CODES = bottle.HTTP_CODES.copy()
HTTP_CODES = dict((y,x) for x,y in HTTP_CODES.iteritems())

app = Bottle()
config = app.config
if not getattr(config, 'load_config'):

        py = sys.version_info
        py3k = py >= (3, 0, 0)

        if py3k:
            from configparser import ConfigParser
        else:
            from ConfigParser import SafeConfigParser as ConfigParser

        from bottle import ConfigDict

        def load_config(self, filename):
            ''' Load values from an *.ini style config file.

                If the config file contains sections, their names are used as
                namespaces for the values within. The two special sections
                ``DEFAULT`` and ``bottle`` refer to the root namespace (no prefix).
            '''
            conf = ConfigParser()
            conf.read(filename)
            for section in conf.sections():
                for key, value in conf.items(section):
                    if section not in ('DEFAULT', 'bottle'):
                        key = section + '.' + key
                    self[key] = value
            return self

        ConfigDict.load_config = load_config

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-c", "--config", dest="config_path",
                  help="path to config file", default="beaker.conf")
    (options, args) = parser.parse_args()
    conf_path = options.config_path
else:
    conf_path = "beaker.conf"

config.load_config(conf_path)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/status')
def status():
    return "This is ChEMBL beaker, version %s" % version

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2smiles/<ctab>')
def ctab2smiles(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return utils._ctab2smiles(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2smiles')
def ctab2smiles():
    data=request.body.getvalue()
    return utils._ctab2smiles(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/smiles2ctab/<smiles>')
def smiles2ctab(smiles):
    data = base64.urlsafe_b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    return utils._smiles2ctab(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/smiles2ctab')
def smiles2ctab():
    data = request.body.getvalue()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    return utils._smiles2ctab(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/inchi2ctab/<inchi>')
def inchi2ctab(inchi):
    inchis = base64.urlsafe_b64decode(inchi)
    return utils._inchi2ctab(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/inchi2ctab')
def inchi2ctab():
    inchis = request.body.getvalue()
    return utils._inchi2ctab(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2inchi/<ctab>')
def ctab2inchi(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return utils._ctab2inchi(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2inchi')
def ctab2inchi():
    data=request.body.getvalue()
    return utils._ctab2inchi(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/inchi2inchiKey/<inchi>')
def inchi2inchiKey(inchi):
    inchis = base64.urlsafe_b64decode(inchi)
    return utils._inchi2inchiKey(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/inchi2inchiKey')
def inchi2inchiKey():
    inchis = request.body.getvalue()
    return utils._inchi2inchiKey(inchis)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2json/<ctab>')
@app.get('/ctab2json/<ctab>/<size>')
@app.get('/ctab2json/<ctab>/<size>/<legend>')
def ctab2json(ctab, size=200, legend=''):
    if not hasattr(utils, '_ctab2json'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'application/json'
    return utils._ctab2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2json')
def ctab2json():
    if not hasattr(utils, '_ctab2json'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.getvalue()
    legend=request.params.get('legend','')
    response.content_type = 'application/json'
    return utils._ctab2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/smiles2json/<ctab>')
@app.get('/smiles2json/<ctab>/<size>')
@app.get('/smiles2json/<ctab>/<size>/<legend>')
def smiles2json(ctab, size=200, legend=''):
    if not hasattr(utils, '_smiles2json'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'application/json'
    return utils._smiles2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/smiles2json')
def smiles2json():
    if not hasattr(utils, '_smiles2json'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.getvalue()
    legend=request.params.get('legend','')
    response.content_type = 'application/json'
    return utils._smiles2json(data, size, legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2svg/<ctab>')
@app.get('/ctab2svg/<ctab>/<size>')
@app.get('/ctab2svg/<ctab>/<size>/<legend>')
def ctab2svg(ctab, size=200, legend=''):
    if not hasattr(utils, '_ctab2svg'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    size = int(size)
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'image/svg+xml'
    return utils._ctab2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2svg')
def ctab2svg():
    if not hasattr(utils, '_ctab2svg'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.getvalue()
    legend=request.params.get('legend','')
    response.content_type = 'image/svg+xml'
    return utils._ctab2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/smiles2svg/<smiles>')
@app.get('/smiles2svg/<smiles>/<size>')
@app.get('/smiles2svg/<smiles>/<size>/<legend>')
def smiles2svg(smiles, size=200, legend=''):
    if not hasattr(utils, '_smiles2svg'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    size = int(size)
    data = base64.urlsafe_b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    response.content_type = 'image/svg+xml'
    return utils._smiles2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/smiles2svg')
def smiles2svg():
    if not hasattr(utils, '_smiles2svg'):
        response.status = HTTP_CODES['Not Implemented']
        return None
    data = request.body.getvalue()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    size = int(request.forms.get('size', 200))
    legend=request.params.get('legend','')
    response.content_type = 'image/svg+xml'
    return utils._smiles2svg(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/ctab2image/<ctab>')
@app.get('/ctab2image/<ctab>/<size>')
@app.get('/ctab2image/<ctab>/<size>/<legend>')
def ctab2image(ctab, size=200, legend=''):
    size = int(size)
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'image/png'
    return utils._ctab2image(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/ctab2image')
def ctab2image():
    size = int(request.forms.get('size', 200))
    data = request.files.values()[0].file.read() if len(request.files) else request.body.getvalue()
    legend=request.params.get('legend','')
    response.content_type = 'image/png'
    return utils._ctab2image(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/smiles2image/<smiles>')
@app.get('/smiles2image/<smiles>/<size>')
@app.get('/smiles2image/<smiles>/<size>/<legend>')
def smiles2image(smiles, size=200, legend=''):
    size = int(size)
    data = base64.urlsafe_b64decode(smiles)
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    response.content_type = 'image/png'
    return utils._smiles2image(data,size,legend)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/smiles2image')
def smiles2image():
    data = request.body.getvalue()
    if not data.startswith('SMILES Name'):
        data = "SMILES Name\n" + data
    size = int(request.forms.get('size', 200))
    legend=request.params.get('legend','')
    response.content_type = 'image/png'
    return utils._smiles2image(data,size,legend)
    
#-----------------------------------------------------------------------------------------------------------------------

@app.get('/image2ctab/<image>')
def image2ctab(image):
    img = base64.urlsafe_b64decode(image)
    return utils._image2ctab(img, config.get('osra_binaries_location', '/usr/bin/osra'))

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/image2ctab')
def image2ctab():
    img = request.body.read()
    return utils._image2ctab(img, config.get('osra_binaries_location', '/usr/bin/osra'))

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/kekulize/<ctab>')
def kekulize(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return utils._kekulize(data)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/kekulize')
def kekulize():
    data = request.body.getvalue()
    return utils._kekulize(data)

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
    return json.dumps(utils._calc(data,lambda x:x.GetNumAtoms()))

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/getNumAtoms')
def getNumAtoms():
    data = request.body.getvalue()
    return json.dumps(utils._calc(data,lambda x:x.GetNumAtoms()))

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/logP/<ctab>')
def logP(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return json.dumps(utils._calc(data, 'MolLogP'))

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/logP')
def logP():
    data = request.body.getvalue()
    return json.dumps(utils._calc(data, 'MolLogP'))

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/TPSA/<ctab>')
def TPSA(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return json.dumps(utils._calc(data, 'TPSA'))

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/TPSA')
def TPSA():
    data = request.body.getvalue()
    return json.dumps(utils._calc(data, 'TPSA'))

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/molWt/<ctab>')
def molWt(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return json.dumps(utils._calc(data, 'MolWt'))
#-----------------------------------------------------------------------------------------------------------------------

@app.post('/molWt')
def molWt():
    data = request.body.getvalue()
    return json.dumps(utils._calc(data, 'MolWt'))

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/descriptors/<ctab>')
def descriptors(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return json.dumps(utils._descriptors(data,request.params))
#-----------------------------------------------------------------------------------------------------------------------

@app.post('/descriptors')
def descriptors():
    data = request.body.getvalue()
    return json.dumps(utils._descriptors(data,request.params))

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/mcs/<ctab>')
def mcs(ctab):
    data = base64.urlsafe_b64decode(ctab)
    return utils._mcs(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/mcs')
def mcs():
    data = request.body.getvalue()
    return utils._mcs(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/smiles2SimilarityMap/<smiles>')
def smiles2SimilarityMap(smiles):
    data = base64.urlsafe_b64decode(smiles)
    response.content_type = 'image/png'
    return utils._smiles2SimilarityMap(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/smiles2SimilarityMap')
def smiles2SimilarityMap():
    data = request.body.getvalue()
    response.content_type = 'image/png'
    return utils._smiles2SimilarityMap(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.get('/sdf2SimilarityMap/<ctab>')
def sdf2SimilarityMap(ctab):
    data = base64.urlsafe_b64decode(ctab)
    response.content_type = 'image/png'
    return utils._sdf2SimilarityMap(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('/sdf2SimilarityMap')
def sdf2SimilarityMap():
    data = request.body.getvalue()
    response.content_type = 'image/png'
    return utils._sdf2SimilarityMap(data,request.params)

#-----------------------------------------------------------------------------------------------------------------------

@app.post('clean')
def clean2D():
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.post('cipStereoInfo')
def stereoInfo():
    pass

#-----------------------------------------------------------------------------------------------------------------------

@app.post('molExport')
def molExport():
    pass

#-----------------------------------------------------------------------------------------------------------------------



if __name__ == "__main__":
    run(app=app, host=config.get('bottle_host', 'localhost'), port=config.get('bottle_port', '8080'),
                                debug=config.get('debug', True), server=config.get('server_middleware', 'tornado'))
else:
    application = app

#-----------------------------------------------------------------------------------------------------------------------
