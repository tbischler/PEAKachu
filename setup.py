try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='PEAKachu',
    version='0.1dev',
    packages=['peakachulib'],
    author='Thorsten Bischler',
    author_email='thorsten.bischler@uni-wuerzburg.de',
    description='Peak calling tool for CLIP-seq data',
    url='',
    install_requires=[
        "biopython >= 1.66",
        "matplotlib >= 1.5.1",
        "pandas >= 0.17.1",
        "pysam >= 0.9.0",
        "bcbio-gff >= 0.6.2",
        "statsmodels >= 0.6.1",
        "numexpr >= 2.5",
        "rpy2 >= 2.7.6"
    ],
    scripts=['bin/peakachu'],
    license='ISC License (ISCL)',
    long_description=open('README.rst').read(),
    classifiers=[
        'License :: OSI Approved :: ISC License (ISCL)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
