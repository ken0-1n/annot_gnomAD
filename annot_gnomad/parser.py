#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .compare import comp_main

def create_parser():
    prog = "annot_gnomad"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
    
    def _create_comp_parser(subparsers):
        
        comp_parser = subparsers.add_parser("comp", help = "compare structural variation results and structural variants in gnomAD")
        comp_parser.add_argument("--input_sv", help = "the result of GenomonSV", type = str, required=True)
        comp_parser.add_argument("--gnomad_vcf", help = "the structural valiants in gnomAD", type = str, required=True)
        comp_parser.add_argument("--out_pref", help = "the output file prefix", type = str, required=True)
        comp_parser.add_argument("--margin", help = "the margin for comparing SVs and SVs in gnomAD", type = int, default = 10)
        
        return comp_parser
        
    comp_parser = _create_comp_parser(subparsers)
    comp_parser.set_defaults(func = comp_main)
    return parser
