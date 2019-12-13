from beaker import config
from beaker.utils import import_class
from beaker.throttle.backends.base import BaseThrottle

throttle = None

throttle_class = config.get('throttle_backend', 'beaker.throttle.backends.cacheThrottle.CacheThrottle')

if throttle_class:
    try:
        throttle = import_class(throttle_class)()
        if not isinstance(throttle, BaseThrottle):
            print("Configured throttle class (%s) is not a BaseThrottle instance, skipping throttling." % throttle_class)
            throttle = None
    except ImportError:
        print('Error importing %s' % throttle_class)






