__author__ = 'mnowotka'

import pickle
import time
from chembl_beaker.beaker.cache.backends.base import BaseCache
from chembl_beaker.beaker.utils.synch import RWLock
from chembl_beaker.beaker import config

#-----------------------------------------------------------------------------------------------------------------------

# Global in-memory store of cache data. Keyed by name, to provide
# multiple named local memory caches.
_caches = {}
_expire_info = {}
_locks = {}

#-----------------------------------------------------------------------------------------------------------------------

class LocMemCache(BaseCache):

    def __init__(self):
        BaseCache.__init__(self)
        name = config.get('cache_loc_mem_instance_name', 'beaker_cache')
        self._cache = _caches.setdefault(name, {})
        self._expire_info = _expire_info.setdefault(name, {})
        self._lock = _locks.setdefault(name, RWLock())

#-----------------------------------------------------------------------------------------------------------------------

    def add(self, key, value):
        key = self.make_key(key)
        pickled = pickle.dumps(value, pickle.HIGHEST_PROTOCOL)
        with self._lock.writer():
            if self._has_expired(key):
                self._set(key, pickled)
                return True
            return False

#-----------------------------------------------------------------------------------------------------------------------

    def get(self, key, default=None):
        key = self.make_key(key)
        pickled = None
        with self._lock.reader():
            if not self._has_expired(key):
                pickled = self._cache[key]
        if pickled is not None:
            try:
                return pickle.loads(pickled)
            except pickle.PickleError:
                return default

        with self._lock.writer():
            try:
                del self._cache[key]
                del self._expire_info[key]
            except KeyError:
                pass
            return default

#-----------------------------------------------------------------------------------------------------------------------

    def _set(self, key, value):
        if len(self._cache) >= self._max_entries:
            self._cull()
        self._cache[key] = value
        self._expire_info[key] = self.get_backend_timeout()

#-----------------------------------------------------------------------------------------------------------------------

    def set(self, key, value):
        key = self.make_key(key)
        pickled = pickle.dumps(value, pickle.HIGHEST_PROTOCOL)
        with self._lock.writer():
            self._set(key, pickled)

#-----------------------------------------------------------------------------------------------------------------------

    def _has_expired(self, key):
        exp = self._expire_info.get(key, -1)
        if exp is None or exp > time.time():
            return False
        return True

#-----------------------------------------------------------------------------------------------------------------------

    def _cull(self):
        if self._cull_frequency == 0:
            self.clear()
        else:
            doomed = [k for (i, k) in enumerate(self._cache) if i % self._cull_frequency == 0]
            for k in doomed:
                self._delete(k)

#-----------------------------------------------------------------------------------------------------------------------

    def _delete(self, key):
        try:
            del self._cache[key]
        except KeyError:
            pass
        try:
            del self._expire_info[key]
        except KeyError:
            pass

#-----------------------------------------------------------------------------------------------------------------------

    def delete(self, key):
        key = self.make_key(key)
        with self._lock.writer():
            self._delete(key)

#-----------------------------------------------------------------------------------------------------------------------

    def clear(self):
        self._cache.clear()
        self._expire_info.clear()

#-----------------------------------------------------------------------------------------------------------------------