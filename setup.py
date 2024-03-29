from setuptools import setup

setup (
    name = "dytop",
    version = "0.1.8",
    author = "Ewerton Rocha Vieira",
    url = "https://github.com/Ewerton-Vieira/dytop.git",
    description = "dytop: combinatorial DYnamics and TOPology",
    long_description = open('README.md').read(),
    ext_package='dytop',
    packages=['dytop'],
    install_requires = ['numpy', 'scipy', 'matplotlib', 'CMGDB', 'pychomp2', 'matplotlib']
)