from setuptools import setup, find_packages

setup(
    name="confWorks",
    version="2025.07",
    description='Python-wrapper for handling conformer ensembles using RDKit and xTB.',
    author='Galymzhan Moldagulov',
    author_email='moldagulovg@gmail.com',
    url='https://github.com/moldagulovg/confWorks',
    packages=['confworks'],
    install_requires=['rdkit', 
                      'xtb', 
                      'tqdm',
                      'numpy',
                      'pandas',
                      'joblib',
                      'chemiscope',
                      'stk',
                      'stko',
                      'opentsne',
                      ]
)
