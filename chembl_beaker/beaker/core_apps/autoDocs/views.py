__author__ = 'mnowotka'

from collections import OrderedDict
from bottle import response, static_file
from chembl_beaker import __version__ as version
import json

from chembl_beaker.beaker import app, config
from chembl_beaker.beaker import STATIC_ROOT
from chembl_beaker.beaker import PARAM_REGEX

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/docs')
def docs():
    return static_file('docs.html', root=STATIC_ROOT)

#-----------------------------------------------------------------------------------------------------------------------

@app.route('/spore', methods = ['OPTIONS', 'GET'])
def spore():
    """Print available functions."""

    ret = {
        "version": version,
        "expected_status": [200],
        "name": "ChEMBL Beaker API live documentation",
        "methods": OrderedDict()
    }

    post_methods = dict()
    get_methods = dict()

    for route in app.routes:
        name = route.name
        method = route.method.upper()
        uname = "%s_%s" % (method, name)
        if not name or method not in ('GET', 'POST'):
            continue
#        if uname in ret["methods"]:
#            #print "powtorka dla %s" % uname
#            params = map(lambda x: x[1:-1].upper(), PARAM_REGEX.findall(route.rule))
#            intersection = list(set(params) & set(ret["methods"][uname]["required_params"]))
#            #print "intersection = %s" % intersection
#            diff = list(set(params) ^ set(ret["methods"][uname]["required_params"]))
#            #print "diff = %s" % diff
#            #print "optional params = %s" % list(set(ret["methods"][uname].get("optional_params", []) + diff))
#            ret["methods"][uname]["required_params"] = intersection
#            ret["methods"][uname]["optional_params"] = list(set(ret["methods"][uname].get("optional_params", []) + diff))
#            continue
        method_info = dict()
        method_info["description"] = route.callback.__doc__
        method_info["method"] = method
        method_info["formats"] = ['text']
        method_info["path"] = PARAM_REGEX.sub(lambda x: ':' + x.group(0)[1:-1].upper(), route.rule)
        method_info["required_params"] = map(lambda x: x[1:-1].upper(), PARAM_REGEX.findall(route.rule))
        if method.upper() == 'POST':
            post_methods[uname] = method_info
        else:
            get_methods[uname] = method_info
    ret["methods"].update([(x, get_methods[x]) for x in sorted(get_methods.keys())])
    ret["methods"].update([(x, post_methods[x]) for x in sorted(post_methods.keys())])
    response.content_type = 'application/json'
    return json.dumps(ret)

#-----------------------------------------------------------------------------------------------------------------------
