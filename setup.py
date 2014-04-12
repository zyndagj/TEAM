#!/usr/bin/python

"""
Setup script for team
"""

from distutils.core import setup, Extension

module1 = Extension('cFetch', sources=["team/cFetch.c"])

setup(name = "team",
	version = "0.1",
	author = "Greg Zynda",
	author_email="gjzynda@indiana.edu",
	license="GNU",
	description = "Indexing of MethylCoder output",
	packages = ["team"],
	ext_modules=[module1])
