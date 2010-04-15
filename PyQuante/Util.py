def parseline(line,format):
    """
    Return the line split and with particular formats applied.

    Format characters:
    x     Skip this record
    s     Convert the record to a string
    f     Convert the record to a float
    d     Convert the record to an int
    i     Convert the record to an int

    Examples:
    >>> parseline('H  0.0 0.0 0.0','sfff')
    'H', 0.0, 0.0, 0.0
    >>> parseline('H  0.0 0.0 0.0','xfff')
    0.0, 0.0, 0.0
    """
    xlat = {'x':None,'s':str,'f':float,'d':int,'i':int}
    result = []
    words = line.split()
    if len(words) < len(format): return None
    for i in xrange(len(format)):
        f = format[i]
        trans = xlat.get(f)
        if trans: result.append(trans(words[i]))
    if len(result) == 0: return None
    if len(result) == 1: return result[0]
    return result

def cleansym(s):
    """This function strips off the garbage (everything after and including
       the first non-letter) in an element name."""
    import re 
    return re.split('[^a-zA-Z]',s)[0]
