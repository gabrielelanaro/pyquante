'''
PyQuante Backend
================

This backend uses another "Handler" class to take care of format file
parsing, these are FormatHandlers.  Each of them is associated with a
format identifier, and each of them parses a string of the given
format with the same naming convention, read and write (They are
really like StringHandlers, but without the format argument)
'''

from FormatHandlers import format_handlers # Defined in __init__.py

from Handlers import StringHandler,FileHandler,FormatUnsupported
import os

class PyQuanteStringHandler(StringHandler):
    def read(self, string, format):
        if not format_handlers.has_key(format):
            raise FormatUnsupported("Format %s not supported"%format)
        fh = format_handlers[format]()
        return fh.read(string)
    def write(self, data, format):
        if not format_handlers.has_key(format):
            raise FormatUnsupported("Format %s not supported"%format)
        fh = format_handlers[format]()
        return fh.write(data)
     
class PyQuanteFileHandler(FileHandler):
    def __init__(self):
        self.string_handler=PyQuanteStringHandler()
    def guess_format(self,filename):
        ext = os.path.splitext(filename)[-1]
        if format_handlers.has_ext(ext):
            return format_handlers.by_ext(ext).key
        else:
            raise FormatUnsupported("Can't recognize the format of this file")

class FormatHandler(object):
    '''
    Handler class used to parse a format.

    Attributes:
    - key: Key used to recognize the format
    - ext: Filename extension used to recognize the format
    - description: Format description, used for documenting
    '''
    key = None
    ext = None
    description = None
    def read(self, string):
        '''
        The functionality is the same as the StringHandler,except for
        the absence of the format argument.
        '''
        return NotImplementedError()
    def write(self, data):
        return NotImplementedError()

