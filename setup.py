import setuptools
import os

def read_version():
    with open(os.path.join('src', 'metaerg', '__init__.py'), 'r', encoding="utf-8") as handle:
        for line in handle:
            if line.startswith('__version__'):
                return line.split('=')[-1].strip().strip('\'')
    raise ValueError("unable to find version")


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='metaerg',
    version=read_version(),
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
    # versions for keras and tensorflow are limited because deepsig is currently not compatible with later versions
    install_requires=['biopython', 'ncbi-datasets-pylib', 'pandas', 'httpx', 'virtualenv', 'h5py', 'pyarrow', 'openpyxl', 'tqdm'], #'keras == 2.15', 'tensorflow == 2.15',
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