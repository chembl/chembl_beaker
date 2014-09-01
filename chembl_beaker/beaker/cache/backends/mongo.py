# -*- coding: utf-8 -*-
# Author Karol Sikora <karol.sikora@laboratorium.ee>, (c) 2012
# Author Michal Nowotka <mmmnow@gmail.com>, (c) 2013-2014

try:
    import cPickle as pickle
except ImportError:
    import pickle
import base64
import pymongo
from datetime import datetime, timedelta
from chembl_beaker.beaker.cache.backends.base import BaseCache
import zlib
import logging
import json
from chembl_beaker.beaker import config

#-----------------------------------------------------------------------------------------------------------------------

MAX_SIZE = 16000000

#-----------------------------------------------------------------------------------------------------------------------

class MongoDBCache(BaseCache):

    def __init__(self):
        super(MongoDBCache, self).__init__()
        self._host = config.get('mongo_host', 'localhost')
        self._port = int(config.get('mongo_port', 27017))
        self._database = config.get('mongo_db', 'beaker_cache')
        self._rshosts = config.get('mongo_rshosts')
        self._rsname = config.get('mongo_rsname')
        self._user = config.get('mongo_user', None)
        self._password = config.get('mongo_pass', None)
        self._socket_timeout_ms = config.get('mongo_socket_simeout_ms', None)
        self._connect_timeout_ms = config.get('mongo_connect_timeout_ms', 300)
        self.compression_level = config.get('mongo_compression_level', 6)
        self._tag_sets = json.loads(config.get('mongo_tag_sets', []))
        self._read_preference = config.get("mongo_read_preference")
        self._collection = config.get('mongo_collection', 'cache')
        self.log = logging.getLogger(__name__)

#-----------------------------------------------------------------------------------------------------------------------

    def add(self, key, value):
        key = self.make_key(key)
        self._base_set('add', key, value)

#-----------------------------------------------------------------------------------------------------------------------

    def set(self, key, value):
        key = self.make_key(key)
        self._base_set('set', key, value)

#-----------------------------------------------------------------------------------------------------------------------

    def _base_set(self, mode, key, value, timeout=None):
        if not timeout:
            timeout = self.default_timeout
        now = datetime.utcnow()
        expires = now + timedelta(seconds=timeout)
        coll = self._get_collection()
        encoded = self._encode(value)
        document_size = len(encoded)
        count = coll.count()
        if count > self._max_entries:
            self._cull()
        data = coll.find_one({'_id': key})
        if data and (mode == 'set' or
                (mode == 'add' and data['expires'] > now)):
            raw = data.get('data')
            if raw and raw == encoded:
                coll.update({'_id': data['_id']}, {'$set': {'expires': expires}}, safe=True)
                return
            else:
                self._delete([key] + data.get('chunks', []))
        if document_size <= MAX_SIZE:
            coll.insert({'_id': key, 'data': encoded, 'expires': expires}, safe=True)
        else:
            chunks = []
            for i in xrange(0, document_size, MAX_SIZE):
                chunk = encoded[i:i+MAX_SIZE]
                aux_key = self.make_key(chunk)
                coll.insert({'_id': aux_key, 'data': chunk}, safe=True)
                chunks.append(aux_key)
            coll.insert({'_id': key, 'chunks': chunks, 'expires': expires}, safe=True)

#-----------------------------------------------------------------------------------------------------------------------

    def _decode(self, data):
        return pickle.loads(zlib.decompress(base64.decodestring(data)))

#-----------------------------------------------------------------------------------------------------------------------

    def _encode(self, data):
        return base64.encodestring(zlib.compress(pickle.dumps(data, pickle.HIGHEST_PROTOCOL), self.compression_level))

#-----------------------------------------------------------------------------------------------------------------------

    def get(self, key, default=None):
        coll = self._get_collection()
        key = self.make_key(key)
        now = datetime.utcnow()
        data = coll.find_one({'_id': key})
        if not data:
            return default
        if data['expires'] < now:
            coll.remove(data['_id'])
            return default
        raw = data.get('data')
        if not raw:
            chunks = data.get('chunks')
            if chunks:
                raw = ''
                for chunk in chunks:
                    raw += coll.find_one({'_id': chunk})['data']
            else:
                return default
        return self._decode(raw)

#-----------------------------------------------------------------------------------------------------------------------

    def delete(self, key):
        key = self.make_key(key)
        coll = self._get_collection()
        data = coll.find_one({'_id': key})
        if data:
            self._delete([key] + data.get('chunks', []))

#-----------------------------------------------------------------------------------------------------------------------

    def _delete(self, ids_to_remove):
        coll = self._get_collection()
        coll.remove({'_id': {'$in':ids_to_remove}})

#-----------------------------------------------------------------------------------------------------------------------

    def clear(self):
        coll = self._get_collection()
        coll.remove({})

#-----------------------------------------------------------------------------------------------------------------------

    def _cull(self):
        if self._cull_frequency == 0:
            self.clear()
            return
        coll = self._get_collection()
        coll.remove({'expires': {'$lte': datetime.utcnow()}})
        #TODO: implement more agressive cull

#-----------------------------------------------------------------------------------------------------------------------

    def _get_collection(self):
        if not getattr(self, '_coll', None):
            self._initialize_collection()
        return self._coll

#-----------------------------------------------------------------------------------------------------------------------

    def _initialize_collection(self):
        try:
            from gevent import monkey
            monkey.patch_socket()
        except ImportError:
            pass

        if self._rsname:
            self.connection = pymongo.MongoReplicaSetClient(self._rshosts, replicaSet=self._rsname,
                read_preference=getattr(pymongo.ReadPreference,self._read_preference, 'PRIMARY'),
                socketTimeoutMS=self._socket_timeout_ms, connectTimeoutMS=self._connect_timeout_ms,
                                                                            tag_sets=self._tag_sets)
        else:
            self.connection = pymongo.Connection(self._host, self._port)

        self._db = self.connection[self._database]
        if self._user and self._password:
            self._db.authenticate(self._user, self._password)
        self._coll= self._db[self._collection]

#-----------------------------------------------------------------------------------------------------------------------