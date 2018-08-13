import sys
import os
from numpy import *

class TinyCache:
    """
    Tiny implementation of data caching
    """

    def __init__(self, cache_path):
        """
        Constructor with set of path to cached data
        """
        self.cache = {}
        self.path = cache_path

    def __contains__(self, key, file_bound=None):
        cache_file = os.path.join(self.path, key)

        if file_bound:
            if os.path.isfile(cache_file) and (os.path.getmtime(cache_file) < os.path.getmtime(cache_file)):
                return True
            else:
                return False
        else:
            if key in self.cache or os.path.isfile(cache_file):
                return True
            else:
                return False

    def get_cache(self, key, file_bound=None):
        cache_file = os.path.join(self.path, key)

        if self.__contains__(key, file_bound):
            if not key in self.cache:
                self.cache[key] = fromfile(cache_file, dtype='float')
            return(self.cache[key])
        else:
            return(None)

    def update_cache(self, key, value):
        cache_file = os.path.join(self.path, key)
        intvalue = asarray(value, dtype='float')
        self.cache[key] = intvalue
        intvalue.tofile(cache_file)
