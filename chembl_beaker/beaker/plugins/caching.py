__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

import base64
import time
import json
from bottle import request, response
from chembl_beaker.beaker import config
from chembl_beaker.beaker.cache import cache

if not cache:
    print "Caching plugin enabled but no cache backend configured, cashing will be skipped..."

if cache and config.get('clear_cache_on_start', False):
    cache.clear()

# ----------------------------------------------------------------------------------------------------------------------


class Caching(object):
    name = 'caching'
    api = 2

    def apply(self, fn, _):
        def _caching(*args, **kwargs):
            start = time.time()
            if not cache:
                from_cache = False
                res = fn(*args, **kwargs)
            else:
                key = json.dumps(args) + json.dumps(kwargs) + json.dumps([(base64.b64encode(k), base64.b64encode(v))
                                                                          for k, v in request.params.items()]) + \
                      request.method + request.path + request.body.read()
                request.body.seek(0)
                try:
                    cached_content, content_type = cache.get(key, (None, None))
                except Exception as e:
                    print e
                    cached_content, content_type = (None, None)
                if cached_content:
                    if content_type:
                        response.headers['Content-Type'] = content_type
                    from_cache = True
                    res = cached_content
                else:
                    from_cache = False
                    res = fn(*args, **kwargs)
                    if type(res) == str:
                        content_type = response.headers.get('Content-Type')
                        try:
                            cache.set(key, (res, content_type))
                        except Exception as e:
                            print e
            if config.get('debug', True):
                end = time.time()
                response.headers['X-ChEMBL-in-cache'] = from_cache
                response.headers['X-ChEMBL-retrieval-time'] = end - start
            return res

        return _caching

# ----------------------------------------------------------------------------------------------------------------------

