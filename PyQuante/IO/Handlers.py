'''
Handlers
========

This module contains a set of superclasses used to handle the reading
and writing of different data from chemical formats.

The implementation of this handlers are in the backends. Currently
there are 2 backends:

- OpenBabelBackend.py 
- PyQuanteBackend.py

These classes are used for example to fetch geometries in Molecule
alternate constructors.
'''


class StringHandler(object):
    """
    StringHandler is an interface for an object that reads data from a
    string and write data in a string
    """
    def read(self, string, format):
        '''
        Reads data from a string written in a given format

        Parameters:
        - string:
        - format: format identifier, a little string

        Raises:
        - FormatUnsupported exception
        '''
        return NotImplementedError()
    def write(self, data, format):
        '''
        Generates a string from the Data object passed, look at read
        for parameters description
        '''
        return NotImplementedError()
    

class FileHandler(object):
    '''
    This class extends the functionality of the StringHandler and uses
    it to read and write the strings from files.
    
    The difference with the StringHandler is the guess_format method,
    used to guess the format from the filename extension.

    Attributes:
    - string_handler: The string Handler that it uses. It's required!!!
    '''
    def read(self, filename, format=None):
        '''
        Read the content of the file "filename".

        Parameters:
        - filename
        - format: if None (the default) attempts to guess the format

        Return:
        - data: a Data instance.
        '''
        if format==None:
            format = self.guess_format(filename)
        string = open(filename,"r").read()
        return self.string_handler.read(string, format)
    def write(self, filename, data, format=None):
        '''
        Write the data (Data instance) in the file specified.
        '''
        if format==None:
            format = self.guess_format(filename)
        fd = open(filename,"w")
        string = self.string_handler.write( data, format)
        fd.write(string)
        fd.close()
    def guess_format(self,filename):
        '''
        Method used to guess the file format from the extension.
        Return:
        - format: the format key if available.
        Raises:
        - FormatUnsupported
        '''
        raise NotImplementedError()

    
class FormatUnsupported(ValueError):
    '''
    Exception raised when something goes wrong with format recognizing
    '''
    pass
