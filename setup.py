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
        "biopython >= 1.74",
        "matplotlib >= 3.1.1",
        "pandas >= 0.25.1",
        "pysam >= 0.15.3",
        "bcbio-gff >= 0.6.6",
        "statsmodels >= 0.10.1",
        "numexpr >= 2.7.0",
        "rpy2 >= 3.1.0"
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
    ext_modules=[Extension("peakachulib.intersection",
                           ["peakachulib/intersection.pyx"])]
)
