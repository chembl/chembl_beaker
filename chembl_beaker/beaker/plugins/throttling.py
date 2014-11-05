__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from datetime import datetime
import pytz
import base64
from bottle import request, response
from chembl_beaker.beaker.throttle import throttle
from chembl_beaker.beaker import config
from chembl_beaker.beaker.utils import import_class

DATE_FORMAT = '%y-%m-%d %H%M%S'

AES = None
secret_key = None
try:
    from Crypto.Cipher import AES
    from Crypto import Random
    secret_key = config.get('throttling_secret_key')
except ImportError:
    pass

#-----------------------------------------------------------------------------------------------------------------------

def verify_key(key):
    if not (AES or secret_key):
        return False
    try:
        iv = Random.new().read(AES.block_size)
        aes = AES.new(secret_key, AES.MODE_CBC, iv)
        now = datetime.utcnow().replace(tzinfo=pytz.utc)
        decoded_key = base64.standard_b64decode(key)
        decrypted_key = aes.decrypt(decoded_key)
        decrypted = filter(bool, decrypted_key.split('\t'))
        dates = [pytz.utc.localize(datetime.strptime(date, DATE_FORMAT)) for date in decrypted]
        if dates[0] < now < dates[1]:
            return True
    except:
        pass

    return False

#-----------------------------------------------------------------------------------------------------------------------

key_verification_function = config.get('throttling_key_verification_function')
if key_verification_function:
    try:
        verified_key = import_class(key_verification_function)
    except ImportError:
        err = 'Error importing function %s' % key_verification_function
        print err
        raise Exception(err)

else:
    if not AES:
        print "API Key verification will be switched off - can't find pycrypto library"
    if not secret_key:
        print "API Key verification will be switched off - no AES key defined"
    verified_key = verify_key

#-----------------------------------------------------------------------------------------------------------------------

def get_identifier(req):
    ip = req.remote_addr
    api_key = req.headers.get('X-ChEMBL-APIKey')
    if api_key:
        if verified_key(api_key):
            return api_key, 'API Key'
        return ip, 'IP (API Key verification failed)'
    return ip, 'IP'

#-----------------------------------------------------------------------------------------------------------------------

class Throttling(object):
    name = 'throttling'
    api = 2

    def apply(self, fn, context):
        def _throttle_check(*args, **kwargs):

            res = fn(*args, **kwargs)
            if not type(res) == str:
                return res
            else:
                identifier, auth_type = get_identifier(request)
                hourly_limit = config.get('throttle_hourly_rate_limit', 60 * 60) if auth_type == 'IP' else \
                                                        config.get('throttle_key_hourly_rate_limit', 60 * 60 * 2)
                daily_limit = config.get('throttle_daily_rate_limit', 60 * 60 * 24) if auth_type == 'IP' else \
                                                    config.get('throttle_key_daily_rate_limit', 60 * 60 * 24 * 2)
                response.headers['X-ChEMBL-Authentication-Type'] = auth_type
                hourly_remaining, daily_remaining = throttle.get_remaining_rates(identifier, auth_type)
                response.headers['X-HourlyRateLimit-Limit'] = hourly_limit
                response.headers['X-DailyRateLimit-Limit'] = daily_limit
                response.headers['X-HourlyRateLimit-Remaining'] = hourly_remaining
                response.headers['X-DailyRateLimit-Remaining'] = daily_remaining
                if not all((hourly_remaining, daily_remaining)):
                    response.status = 429
                    response.body = 'Too many requests, try again later'
                    return response

                throttle.accessed(identifier, auth_type)

                return res

        return _throttle_check

#-----------------------------------------------------------------------------------------------------------------------