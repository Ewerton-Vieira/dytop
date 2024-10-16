from setuptools import setup, find_packages

setup (
    name = "dytop",
    version = "0.1.13",
    author = "Ewerton Rocha Vieira",
    url = "https://github.com/Ewerton-Vieira/dytop.git",
    description = "dytop: combinatorial DYnamics and TOPology",
    long_description = open('README.md').read(),
    long_description_content_type='text/markdown',
    ext_package='dytop',
    packages=find_packages(),
    install_requires = ['numpy', 'scipy', 'matplotlib', 'CMGDB', 'pychomp2']
)