#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Global settings and constants.

Directories can be overridden via the environment variables
``SPLICE_GTF_CACHE_DIR`` / ``SPLICE_ALLGTF_DIR``.
"""

import os

GTF_CACHE_DIR = os.environ.get('SPLICE_GTF_CACHE_DIR', './gtf_cache')
ALLGTF_DIR = os.environ.get('SPLICE_ALLGTF_DIR', '/mnt/g/bioinfo/refer/Allgtf')

# Maximum 3'UTR extension downstream of the CDS when searching for a stop codon
ALT_DOWNSTREAM_EXTENSION_BP = 1500

# Standard rMATS output file names
RMATS_EVENT_FILES = {
    'SE': 'SE.MATS.JC.txt',
    'RI': 'RI.MATS.JC.txt',
    'A3SS': 'A3SS.MATS.JC.txt',
    'A5SS': 'A5SS.MATS.JC.txt',
    'MXE': 'MXE.MATS.JC.txt',
}


def get_cache_dir() -> str:
    """Return the GTF cache directory, creating it if necessary."""
    os.makedirs(GTF_CACHE_DIR, exist_ok=True)
    return GTF_CACHE_DIR
