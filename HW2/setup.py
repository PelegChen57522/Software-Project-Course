from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['kmeansmodule.c'])
setup(name='mykmeanssp',
     version='1.0',
     description='K-means algo - HW2',
     ext_modules=[module])
