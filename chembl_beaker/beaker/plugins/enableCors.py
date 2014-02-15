__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from bottle import request, response
from chembl_beaker.beaker import config

#-----------------------------------------------------------------------------------------------------------------------

class EnableCors(object):
    name = 'enable_cors'
    api = 2

    def apply(self, fn, context):
        def _enable_cors(*args, **kwargs):
            # set CORS headers
            response.headers['Access-Control-Allow-Origin'] = config.get('access_control_allow_origin', '*')
            response.headers['Access-Control-Allow-Methods'] = config.get('access_control_allow_methods',
                                                                                        'GET, POST, PUT, OPTIONS')
            response.headers['Access-Control-Allow-Headers'] = config.get('access_control_allow_headers',
                                                    'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token')

            if request.method != 'OPTIONS':
                # actual request; reply with the actual response
                return fn(*args, **kwargs)

        return _enable_cors

#-----------------------------------------------------------------------------------------------------------------------