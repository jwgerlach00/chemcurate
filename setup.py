from setuptools import setup, find_packages


setup(
    name='chemcurate',
    version='0.0.1',
    license='MIT',
    author='Jacob Gerlach',
    author_email='jwgerlach00@gmail.com',
    url='https://github.com/jwgerlach00/chemcurate',
    description='Autocuration of small-molecule data from databases.',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'assets': [
        'uniprot_mapping.yaml'
    ]},
    python_requires='>=3.7',
    install_requires=[
        'pandas',
        'PyYAML',
        'requests',
        'setuptools'
    ],
)