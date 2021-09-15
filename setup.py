from setuptools import setup, Extension
setup(
    name='pyRSD_CoEv',
    url='',
    description=('An python package for selection sweep detection and co-evolutionary gene cluster identification.'),
    #ext_modules=cythonize(ext_modules),    
    package_dir={'': 'pyRSD_CoEv'},
    packages=['pyRSD_CoEv'],
    install_requires=[
        'numpy',
        'scipy',
        'networkx',
        #'multiprocessing',
    ],
    dependency_links=[
        ''
    ],
    scripts=['bin/convertFormat','bin/caculateRsd','bin/mergeANDsmooth','bin/annotationRSD','bin/coevCluster','bin/plot']
)