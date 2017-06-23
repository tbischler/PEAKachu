try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from distutils.extension import Extension

setup(
    name='PEAKachu',
    version='0.1.0',
    packages=['peakachulib'],
    author='Thorsten Bischler',
    author_email='thorsten.bischler@uni-wuerzburg.de',
    description='Peak calling tool for CLIP-seq data',
    url='',
    install_requires=[
        "biopython >= 1.69",
        "matplotlib >= 2.0.2",
        "pandas >= 0.20.2",
        "pysam >= 0.11.2.2",
        "bcbio-gff >= 0.6.4",
        "statsmodels >= 0.8.0",
        "numexpr >= 2.6.2",
        "rpy2 >= 2.8.5"
    ],
    scripts=['bin/peakachu'],
    license='ISC License (ISCL)',
    long_description=open('README.rst').read(),
    classifiers=[
        'License :: OSI Approved :: ISC License (ISCL)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules = [Extension("peakachulib.intersection", [ "peakachulib/intersection.pyx" ])]
)
