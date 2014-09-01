"Memcached cache backend"

#-----------------------------------------------------------------------------------------------------------------------

import json
try:
    import cPickle as pickle
except ImportError:
    import pickle

from chembl_beaker.beaker.cache.backends.base import BaseCache
from chembl_beaker.beaker import config
from chembl_beaker.beaker.utils.functional import cached_property

#-----------------------------------------------------------------------------------------------------------------------

class BaseMemcachedCache(BaseCache):
    def __init__(self, library, value_not_found_exception):
        super(BaseMemcachedCache, self).__init__()
        if config.get('memcached_servers'):
            self._servers = json.loads(config.get('memcached_servers'))
        else:
            self._servers = [config.get('memcached_server', '127.0.0.1:11211')]

        # The exception type to catch from the underlying library for a key
        # that was not found. This is a ValueError for python-memcache,
        # pylibmc.NotFound for pylibmc, and cmemcache will return None without
        # raising an exception.
        self.LibraryValueNotFoundException = value_not_found_exception

        self._lib = library

#-----------------------------------------------------------------------------------------------------------------------

    @property
    def _cache(self):
        """
        Implements transparent thread-safe access to a memcached client.
        """
        if getattr(self, '_client', None) is None:
            self._client = self._lib.Client(self._servers)

        return self._client

#-----------------------------------------------------------------------------------------------------------------------

    def add(self, key, value):
        key = self.make_key(key)
        return self._cache.add(key, value, self.default_timeout)

#-----------------------------------------------------------------------------------------------------------------------

    def get(self, key, default=None):
        key = self.make_key(key)
        val = self._cache.get(key)
        if val is None:
            return default
        return val

#-----------------------------------------------------------------------------------------------------------------------

    def set(self, key, value):
        key = self.make_key(key)
        self._cache.set(key, value, self.default_timeout)

#-----------------------------------------------------------------------------------------------------------------------

    def delete(self, key):
        key = self.make_key(key)
        self._cache.delete(key)

#-----------------------------------------------------------------------------------------------------------------------

    def close(self, **kwargs):
        self._cache.disconnect_all()

#-----------------------------------------------------------------------------------------------------------------------

    def clear(self):
        self._cache.flush_all()


#-----------------------------------------------------------------------------------------------------------------------

class MemcachedCache(BaseMemcachedCache):
    "An implementation of a cache binding using python-memcached"
    def __init__(self):
        import memcache
        super(MemcachedCache, self).__init__(library=memcache, value_not_found_exception=ValueError)

    @property
    def _cache(self):
        if getattr(self, '_client', None) is None:
            self._client = self._lib.Client(self._servers, pickleProtocol=pickle.HIGHEST_PROTOCOL)
        return self._client


#-----------------------------------------------------------------------------------------------------------------------

class PyLibMCCache(BaseMemcachedCache):
    "An implementation of a cache binding using pylibmc"
    def __init__(self):
        import pylibmc
        super(PyLibMCCache, self).__init__(library=pylibmc, value_not_found_exception=pylibmc.NotFound)

    @cached_property
    def _cache(self):
        client = self._lib.Client(self._servers)
        return client

#-----------------------------------------------------------------------------------------------------------------------
