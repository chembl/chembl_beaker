__author__ = 'mnowotka'

from chembl_beaker.beaker import config

#-----------------------------------------------------------------------------------------------------------------------

class BaseThrottle(object):

    def __init__(self):
        self.hourly_rate_limit = config.get('throttle_hourly_rate_limit', 60 * 60)
        self.daily_rate_limit = config.get('throttle_daily_rate_limit', 60 * 60 * 24)
        self.key_hourly_rate_limit = config.get('throttle_key_hourly_rate_limit', 60 * 60 * 2)
        self.key_daily_rate_limit = config.get('throttle_key_daily_rate_limit', 60 * 60 * 24 * 2)

    def get_remaining_rates(self, identifier, type='IP'):
        raise NotImplementedError('subclasses of BaseThrottle must provide a get_remaining_rates method')

    def accessed(self, identifier, type='IP'):
        raise NotImplementedError('subclasses of BaseThrottle must provide an accessed method')

#-----------------------------------------------------------------------------------------------------------------------
