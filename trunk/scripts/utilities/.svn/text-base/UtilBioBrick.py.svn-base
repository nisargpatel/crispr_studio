#!/usr/bin/env python

"""
Implements basic functionality for dealing with BioBrick assembly.
"""

__author__ = "Ruben Acuna"
__copyright_ = "Copyright(c) 2011, ASU iGEM Team"

################################## CONSTANTS ###################################
BIOBRICK_SUFFIX = "tactagtagcggccgctgcag"
BIOBRICK_PREFIX = "gaattcgcggccgcttctagag"
BIOBRICK_PREFIX_ATG = "gaattcgcggccgcttctag"

################################## FUNCTIONS ###################################
def makeBioBrick(sequence):
    """Given a sequence, concatenates the appropriate BioBrick prefix and suffix."""

    if sequence.startswith("atg"):
        return BIOBRICK_PREFIX_ATG + sequence + BIOBRICK_SUFFIX
    else:
        return BIOBRICK_PREFIX + sequence + BIOBRICK_SUFFIX