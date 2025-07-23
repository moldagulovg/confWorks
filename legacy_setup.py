from setuptools import setup, find_packages

setup(
    name="confWorks",
    version="2025.07",
    description='Convenient python-wrapper for handling multconformer RDKit molecules and conducting geometry optimization and conformational sampling routines (GFN-xTB semi-empirical theory levels).',
    author='Galymzhan Moldagulov',
    author_email='moldagulovg@gmail.com',
    url='https://github.com/moldagulovg/confWorks',
    packages=['confworks'],
    install_requires=['rdkit==2025.3.3',
                      'tqdm',
                      'numpy==1.26.4',
                      'pandas',
                      'joblib',
                      'chemiscope',
                      'stk',
                      'stko',
                      'opentsne',
                      'py3Dmol',
                      ]
)
