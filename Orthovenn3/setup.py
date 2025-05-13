from setuptools import setup, find_packages

setup(
    name='Orthovenn3_Orthogroup_Filter',
    version='0.2.0',
    author='Gabriel Dall\'Alba',
    author_email='gdalba.biol@gmail.com',
    description='A project to filter entries from an Excel file based on identifiers from a text file. Designed for Orthovenn3 output.',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'pandas',
        'openpyxl',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'orthovenn3_orthogroup_filter=main:main',
        ],
    },
)