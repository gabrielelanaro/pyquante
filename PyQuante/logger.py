'''
This module sets the output properties and configurations
'''
import logging
import sys
# Output logger creation and initialization
output = logging.getLogger("pyquante")
output.propagate = 0 # Avoid propagator to the root logger, may interfere with other applications

class NullHandler(logging.Handler):
    '''
    Handler that does nothing
    '''
    def emit(self,record):
        pass
# adding a NullHandler to prevent errors like "hanler not found"
output.addHandler(NullHandler())

configured_handlers = []

def configure_output(filename=None, stream=sys.stdout, level=logging.INFO):
    """
    Shortcut function that configures output handling of pyquante
    output with reasonable defaults
    
    Arguments:
    - filename:
    """
    # Delete previous configurations
    global configured_handlers
    for handler in configured_handlers :
        output.removeHandler(handler)
    configured_handlers = []

    # Default formatter, just the message
    default_formatter = logging.Formatter("%(message)s")
    
    # Process options
    if filename:
        h = logging.FileHandler(filename)
        h.setFormatter(default_formatter)
        output.addHandler(h)
        configured_handlers.append(h)
    if stream:
        h = logging.StreamHandler(stream)
        h.setFormatter(default_formatter)
        output.addHandler(h)
        configured_handlers.append(h)
    output.setLevel(level)
