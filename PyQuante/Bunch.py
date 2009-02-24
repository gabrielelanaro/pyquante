"""\
 Bunch.py: Generic collection that can be accessed as a class. This can
  be easily overrided and used as container for a variety of objects.

 This program is part of the PyQuante quantum chemistry program suite

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Bunch is a generic object that can hold data. You can access it
#  using object syntax: b.attr1, b.attr2, etc., which makes it simpler
#  than a dictionary

class Bunch:
    def __init__(self,**keyw):
        for key in keyw.keys():
            self.__dict__[key] = keyw[key]
        return
