from setuptools import Extension, setup

setup(
    name="symnmfModule",
    version="1.0",
    description="Final Project",
    ext_modules=[
        Extension("symnmfModule", sources=[ "symnmf.c" , "symnmfmodule.c"]),
    ],
)

