try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='PEAKachu',
    version='0.1dev',
    packages=['peakachulib'],
    author='Thorsten Bischler',
    author_email='',
    description='Peak calling tool for CLIP-seq data',
    url='',
    install_requires=[
        "pysam >= 0.8.3"
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
