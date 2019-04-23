__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import json
from bottle import BaseRequest, request, response
from netaddr import IPNetwork, IPAddress
from chembl_beaker.beaker import config

#-----------------------------------------------------------------------------------------------------------------------

IP_WHITELIST = config.get('ip_whitelist')

try:
    if IP_WHITELIST:
        IP_WHITELIST = [IPNetwork(ip) for ip in json.loads(IP_WHITELIST)]
except:
    pass

IP_BLACKLIST = config.get('ip_blacklist')

try:
    if IP_BLACKLIST:
        IP_BLACKLIST = [IPNetwork(ip) for ip in json.loads(IP_BLACKLIST)]
except:
    pass

#-----------------------------------------------------------------------------------------------------------------------

class Restrictions(object):
    name = 'restrictions'
    api = 2

    def __init__(self):
        if config.get('request_max_size'):
            BaseRequest.MEMFILE_MAX = int(config['request_max_size'])

    def apply(self, fn, _):
        def _check_restrictions(*args, **kwargs):
            if config.get('request_max_size') and request.content_length > int(config['request_max_size']):
                response.status = 413
                response.body = 'Request size larger than %s bytes' % config['request_max_size']
                return response

            if request.remote_addr:
                ip = IPAddress(request.remote_addr)
                if IP_WHITELIST:
                    if not any([ip in white for white in IP_WHITELIST]):
                        response.status = 403
                        return response
                elif IP_BLACKLIST:
                    if any([ip in black for black in IP_BLACKLIST]):
                        response.status = 403
                        return response

            return fn(*args, **kwargs)

        return _check_restrictions

#-----------------------------------------------------------------------------------------------------------------------