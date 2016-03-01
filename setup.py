#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import os

install_requires = ['numpy>=1.3.0',]

setup(
    name = "neuropype_graph",
    version = '0.0.1dev',
    packages = ['neuropype_graph'],
    install_requires=install_requires,
    author = "David Meunier",
    description = "Graph analysis for neuropype (using Nipype, and neuropype_ephy)"
)

