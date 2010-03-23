'''
This module contains a class that is useful to the PyQuanteBackend
'''
class HandlerList(object):
    def __init__(self,*list):
        self.list = list
    def has_key(self,key):
        for hand in self.list:
            if hand.key == key:
                return True
        return False
    def has_ext(self,ext):
        for hand in self.list:
            if hand.ext == ext:
                return True
        return False
    def by_ext(self,ext):
        for hand in self.list:
            if hand.ext == ext:
                return hand
        raise KeyError()
    def __getitem__(self,key):
        for hand in self.list:
            if hand.key == key:
                return hand
        raise KeyError(key)
