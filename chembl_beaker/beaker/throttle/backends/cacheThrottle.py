__author__ = 'mnowotka'

import pytz
from datetime import datetime, timedelta
from chembl_beaker.beaker.throttle.backends.base import BaseThrottle
from chembl_beaker.beaker.cache import cache

#-----------------------------------------------------------------------------------------------------------------------

class CacheThrottle(BaseThrottle):

    def __init__(self):
        BaseThrottle.__init__(self)
        if not cache:
            err = "CacheThrottle class can't work without cache..."
            print err
            raise Exception(err)

    def get_remaining_rates(self, identifier, type='IP'):
        start_usage = (None, self.hourly_rate_limit, self.daily_rate_limit) if type == 'IP' else \
                                                        (None, self.key_hourly_rate_limit, self.key_daily_rate_limit)
        try:                                                
            usage = cache.get(identifier, start_usage)
        except:
            usage = start_usage
        return usage[1], usage[2]

    def accessed(self, identifier, type='IP'):
        """
        Handles recording the user's access.

        Stores the current timestamp in the "accesses" list within the cache.
        """
        start_usage = (None, self.hourly_rate_limit, self.daily_rate_limit) if type == 'IP' else \
                                                        (None, self.key_hourly_rate_limit, self.key_daily_rate_limit)
        _, hourly_rate_limit, daily_rate_limit = start_usage
        try:
            usage = cache.get(identifier, start_usage)
        except:
            usage = start_usage
        now = datetime.utcnow().replace(tzinfo=pytz.utc)
        last_usage, old_hourly_rate_limit, old_daily_rate_limit = usage
        if not last_usage:
            try:
                cache.set(identifier, (now, max(old_hourly_rate_limit-1, 0), max(old_daily_rate_limit-1, 0)))
            except:
                pass
        else:
            if last_usage + timedelta(hours=1) >= now and now.hour == last_usage.hour:
                new_hourly_rate_limit = max(old_hourly_rate_limit-1, 0)
            else:
                new_hourly_rate_limit = hourly_rate_limit
            if last_usage + timedelta(days=1) >= now and now.day == last_usage.day:
                new_daily_rate_limit = max(old_daily_rate_limit-1, 0)
            else:
                new_daily_rate_limit = daily_rate_limit
            try:    
                cache.set(identifier, (now, new_hourly_rate_limit, new_daily_rate_limit))
            except:
                pass

#-----------------------------------------------------------------------------------------------------------------------
