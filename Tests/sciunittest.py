#!/usr/bin/env python
import unittest

class TestCase(unittest.TestCase):
    def assertInside(self, first, second, error, msg=None):
        """Fail is the second number isn't within a certain error of the first."""
        if not (second-error) < first < (second+error):
            raise self.failureException, (msg or '%r != %r (+-%r)' % (first, second, error))
