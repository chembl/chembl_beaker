__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

import time
import hashlib
from chembl_beaker.beaker.utils import import_class
from chembl_beaker.beaker import config

#-----------------------------------------------------------------------------------------------------------------------

def default_key_func(key, key_prefix):
    """
    Default function to generate keys.

    Constructs the key used by all other methods. By default it prepends
    the `key_prefix'. KEY_FUNCTION can be used to specify an alternate
    function with custom key making behavior.
    """
    return hashlib.md5(':'.join([key_prefix, key])).hexdigest()

#-----------------------------------------------------------------------------------------------------------------------

def get_key_func(key_func):
    """
    Function to decide which key function to use.

    Defaults to ``default_key_func``.
    """
    if key_func is not None:
        if callable(key_func):
            return key_func
        else:
            return import_class(key_func)
    return default_key_func

#-----------------------------------------------------------------------------------------------------------------------

class BaseCache(object):

    def __init__(self):
        timeout = config.get('cache_timeout', 300)
        if timeout is not None:
            try:
                timeout = int(timeout)
            except (ValueError, TypeError):
                timeout = 300
            self.default_timeout = timeout
        max_entries = config.get('cache_max_entries', 300000)
        try:
            self._max_entries = int(max_entries)
        except (ValueError, TypeError):
            self._max_entries = 300
        cull_frequency = config.get('cache_cull_frequency', 3)
        try:
            self._cull_frequency = int(cull_frequency)
        except (ValueError, TypeError):
            self._cull_frequency = 3
        self.key_prefix = config.get('cache_key_prefix', '')
        self.key_func = get_key_func(config.get('cache_key_function', None))

#-----------------------------------------------------------------------------------------------------------------------

    def make_key(self, key):
        """Constructs the key used by all other methods. By default it
        uses the key_func to generate a key (which, by default,
        prepends the `key_prefix' and 'version'). An different key
        function can be provided at the time of cache construction;
        alternatively, you can subclass the cache backend to provide
        custom key making behavior.
        """

        new_key = self.key_func(key, self.key_prefix)
        return new_key

#-----------------------------------------------------------------------------------------------------------------------

    def add(self, key, value):
        """
        Set a value in the cache if the key does not already exist. If
        timeout is given, that timeout will be used for the key; otherwise
        the default cache timeout will be used.
        Returns True if the value was stored, False otherwise.
        """
        raise NotImplementedError('subclasses of BaseCache must provide an add() method')

#-----------------------------------------------------------------------------------------------------------------------

    def get(self, key, default=None):
        """
        Fetch a given key from the cache. If the key does not exist, return
        default, which itself defaults to None.
        """
        raise NotImplementedError('subclasses of BaseCache must provide a get() method')

#-----------------------------------------------------------------------------------------------------------------------

    def set(self, key, value):
        """
        Set a value in the cache. If timeout is given, that timeout will be
        used for the key; otherwise the default cache timeout will be used.
        """
        raise NotImplementedError('subclasses of BaseCache must provide a set() method')

#-----------------------------------------------------------------------------------------------------------------------

    def delete(self, key):
        """
        Delete a key from the cache, failing silently.
        """
        raise NotImplementedError('subclasses of BaseCache must provide a delete() method')

#-----------------------------------------------------------------------------------------------------------------------

    def clear(self):
        """Remove *all* values from the cache at once."""
        raise NotImplementedError('subclasses of BaseCache must provide a clear() method')

#-----------------------------------------------------------------------------------------------------------------------

    def get_backend_timeout(self):
        """
        Returns the timeout value usable by this backend based upon the provided
        timeout.
        """
        return time.time() + self.default_timeout

#-----------------------------------------------------------------------------------------------------------------------