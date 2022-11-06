import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='metaerg',
    version='2.2.32',
    packages=setuptools.find_packages(where='src'),
    url='https://github.com/kinestetika/MetaErg',
    license='MIT',
    author='Marc Strous',
    author_email='mstrous@ucalgary.ca',
    description='Annotation of genomes and contigs',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Natural Language :: English',
                 'Operating System :: OS Independent',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.10',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='repeat-regions genes functions taxonomy',
    project_urls={'Source': 'https://github.com/kinestetika/MetaErg'},
    package_dir={'': 'src'},
    package_data={'': ['functional_gene_data']},
    python_requires='>=3.10',
    install_requires=['biopython', 'ncbi-datasets-pylib', 'pandas', 'httpx', 'virtualenv', 'h5py', 'pyarrow'],
    extras_require={  # Optional
        'dev': ['setuptools', 'build', 'twine'],
        'test': []
    },
    entry_points={  # Optional
        'console_scripts': [
            'metaerg=metaerg.main:main',
        ],
    }
)