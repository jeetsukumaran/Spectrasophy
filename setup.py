#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(
    name="spectrasophy",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=["spectrasophy"],
    scripts=[
        "bin/spectrasophy-simulate.py",
        "bin/spectrasophy-normalize.py",
        "bin/spectrasophy-reject.py",
        ],
    url="http://pypi.python.org/pypi/spectrasophy/",
    test_suite = "spectrasophy.test",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.txt").read(),
    # install_requires=[ ],
)
