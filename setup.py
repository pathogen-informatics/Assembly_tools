import os
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='Assembly_tools',
    version='0.1.1',
    description='Scripts relating to genome assembly',
    long_description=read('README.md'),
    packages = find_packages(),
    author='Martin Hunt',
    author_email='mh12@sanger.ac.uk',
    url='https://github.com/martinghunt/Assembly_tools',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=['nose >= 1.3','pyfastaq >= 3.2.0','pysam >= 0.7.5'],
    license='GPLv3',
)
