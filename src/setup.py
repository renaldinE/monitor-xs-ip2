from setuptools import setup, find_packages

setup(
    name='monitor_xs_ip2',
    version='1.0.0',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'pandas',
        'scipy',
        'numpy',
        'matplotlib',
        'xlsxwriter',
        'periodictable',
        'math',
        'pysrim'
    ],
)
