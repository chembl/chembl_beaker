__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

from datetime import datetime
import pytz
import base64
from chembl_beaker.beaker import config

DATE_FORMAT = '%y-%m-%d %H%M%S'

AES = None
try:
    from Crypto.Cipher import AES
    secret_key = config.get('throttling_secret_key')
except ImportError:
    pass

#-----------------------------------------------------------------------------------------------------------------------

def generate_key(validity):
    if secret_key:
        aes = AES.new(secret_key, AES.MODE_CBC)
    now = datetime.utcnow().replace(tzinfo=pytz.utc)
    valid_from = now
    valid_to = now + validity
    input_string = "%s\t\t%s" % tuple([date.strftime(DATE_FORMAT) for date in (valid_from, valid_to)])
    return base64.standard_b64encode(aes.encrypt(input_string))

#-----------------------------------------------------------------------------------------------------------------------
