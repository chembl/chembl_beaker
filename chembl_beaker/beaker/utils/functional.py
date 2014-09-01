__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

def _call(mols, fn, *args, **kwargs):
    return [getattr(m, fn)(*args, **kwargs) for m in mols if m and hasattr(m, fn) ]

#-----------------------------------------------------------------------------------------------------------------------

def _apply(mols, fn, *args, **kwargs):
    return [fn(m, *args, **kwargs) for m in mols if m ]

#-----------------------------------------------------------------------------------------------------------------------

class cached_property(object):
    """
    Decorator that converts a method with a single self argument into a
    property cached on the instance.

    Optional ``name`` argument allows you to make cached properties of other
    methods. (e.g.  url = cached_property(get_absolute_url, name='url') )
    """
    def __init__(self, func, name=None):
        self.func = func
        self.name = name or func.__name__

    def __get__(self, instance, type=None):
        if instance is None:
            return self
        res = instance.__dict__[self.name] = self.func(instance)
        return res

#-----------------------------------------------------------------------------------------------------------------------