__author__ = 'mnowotka'

#-----------------------------------------------------------------------------------------------------------------------

try:
    import cPickle as pickle
except ImportError:
    import pickle
from threading import Thread
from Queue import Queue
import base64
import sqlite3
from datetime import datetime
from chembl_beaker.beaker.cache.backends.base import BaseCache
from chembl_beaker.beaker import config
import pytz

#-----------------------------------------------------------------------------------------------------------------------

def value_to_db_datetime(value):
    if value is None:
        return None
    if value.tzinfo is not None and value.tzinfo.utcoffset(value) is not None:
        value = value.astimezone(pytz.utc).replace(tzinfo=None)
    return value

#-----------------------------------------------------------------------------------------------------------------------

class MultiThreadOK(Thread):
    def __init__(self, db):
        super(MultiThreadOK, self).__init__()
        self.db=db
        self.reqs=Queue()
        self.start()
    def run(self):
        cnx = sqlite3.connect(self.db)
        cursor = cnx.cursor()
        while True:
            req, arg, res = self.reqs.get()
            if req=='--close--': break
            cursor.execute(req, arg)
            if res:
                for rec in cursor:
                    res.put(rec)
                res.put('--no more--')
        cnx.close()
    def execute(self, req, arg=None, res=None):
        self.reqs.put((req, arg or tuple(), res))
    def select(self, req, arg=None):
        res=Queue()
        self.execute(req, arg, res)
        while True:
            rec=res.get()
            if rec=='--no more--': break
            yield rec
    def close(self):
        self.execute('--close--')

#-----------------------------------------------------------------------------------------------------------------------

class SQLiteCache(BaseCache):
    def __init__(self):
        super(SQLiteCache, self).__init__()
        self._table = config.get('sqlite_cache_table', 'beaker_cache')
        self.sql = MultiThreadOK(config.get('sqlite_dbfile', ':memory:'))
        self.sql.execute("CREATE TABLE IF NOT EXISTS %s (cache_key varchar(255) PRIMARY KEY NOT NULL, value text, expires datetime NOT NULL)" % self._table)

#-----------------------------------------------------------------------------------------------------------------------

    def get(self, key, default=None):
        key = self.make_key(key)

        try:
            row = next(self.sql.select("SELECT cache_key, value, expires FROM %s "
                               "WHERE cache_key = ?" % self._table, [key]))
        except StopIteration:
            row = None

        if row is None:
            return default
        now = datetime.utcnow().replace(tzinfo=pytz.utc)
        expires = pytz.utc.localize(datetime.strptime(row[2], '%Y-%m-%d %H:%M:%S'))

        if expires < now:
            self.sql.execute("DELETE FROM %s "
                               "WHERE cache_key = ?" % self._table, [key])
            return default
        value = row[1]
        return pickle.loads(base64.b64decode(value))

#-----------------------------------------------------------------------------------------------------------------------

    def set(self, key, value):
        key = self.make_key(key)
        self._base_set('set', key, value, self.get_backend_timeout())

#-----------------------------------------------------------------------------------------------------------------------

    def add(self, key, value):
        key = self.make_key(key)
        return self._base_set('add', key, value, self.get_backend_timeout())

#-----------------------------------------------------------------------------------------------------------------------

    def _base_set(self, mode, key, value, timeout):

        row = next(self.sql.select("SELECT COUNT(*) FROM %s" % self._table))
        num = row[0]
        now = datetime.utcnow().replace(tzinfo=pytz.utc)
        now = now.replace(microsecond=0)
        if timeout is None:
            exp = datetime.max
        else:
            exp = datetime.fromtimestamp(timeout)
        exp = exp.replace(microsecond=0)
        if num > self._max_entries:
            self._cull(now)
        pickled = pickle.dumps(value, pickle.HIGHEST_PROTOCOL)
        b64encoded = base64.b64encode(pickled)

        # Note: typecasting for datetimes is needed by some 3rd party
        # database backends. All core backends work without typecasting,
        # so be careful about changes here - test suite will NOT pick
        # regressions.

        try:
            result = next(self.sql.select("SELECT cache_key, expires FROM %s WHERE cache_key = ?" % self._table, [key]))
        except StopIteration:
            result = None
        if result:
            current_expires = result[1]

        exp = value_to_db_datetime(exp)
        if result and (mode == 'set' or (mode == 'add' and current_expires < now)):
            self.sql.execute("UPDATE %s SET value = ?, expires = ? WHERE cache_key = ?" % self._table,
                           (b64encoded, exp, key))
        else:
            self.sql.execute("INSERT INTO %s (cache_key, value, expires) VALUES (?, ?, ?)" % self._table,
                           (key, b64encoded, exp))

            return True

#-----------------------------------------------------------------------------------------------------------------------

    def delete(self, key):
        key = self.make_key(key)
        self.sql.execute("DELETE FROM %s WHERE cache_key = ?" % self._table, [key])

#-----------------------------------------------------------------------------------------------------------------------

    def _cull(self, now):
        if self._cull_frequency == 0:
            self.clear()
        else:
            # When USE_TZ is True, 'now' will be an aware datetime in UTC.
            now = now.replace(tzinfo=None)
            self.sql.execute("DELETE FROM %s WHERE expires < ?" % self._table,
                           [value_to_db_datetime(now)])
            row = next(self.sql.select("SELECT COUNT(*) FROM %s" % self._table))
            num = row[0]
            if num > self._max_entries:
                cull_num = num // self._cull_frequency
                row = next(self.sql.select(
                    "SELECT cache_key FROM %s ORDER BY cache_key LIMIT 1 OFFSET %s" % (self._table,
                    cull_num)))
                self.sql.execute("DELETE FROM %s "
                               "WHERE cache_key < ?" % self._table,
                               [row[0]])

#-----------------------------------------------------------------------------------------------------------------------

    def clear(self):
        self.sql.execute('DELETE FROM %s' % self._table)

#-----------------------------------------------------------------------------------------------------------------------