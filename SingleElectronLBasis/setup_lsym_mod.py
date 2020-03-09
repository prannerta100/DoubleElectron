from distutils.core import setup, Extension
setup(name = 'LSymBasisModule', version = '1.0',  \
   ext_modules = [Extension('LSymBasisModule', ['lsym_mod.c'])])
