from setuptools import setup, find_packages

setup (
    name = "dytop",
    version = "0.1",
    # packages = find_packages(),
    author = "Ewerton Rocha Vieira",
    url = "https://github.com/Ewerton-Vieira/dytop.git",
    description = "dytop: combinatorial DYnamics and TOPology",
    long_description = open('README.md').read(),
    package_dir = {'': 'src'},
    ext_package='dytop',
    packages=['dytop'],
    install_requires = ['numpy', 'scipy', 'matplotlib', 'CMGDB'],
)