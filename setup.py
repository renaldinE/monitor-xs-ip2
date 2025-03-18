from setuptools import setup, find_packages

setup(
    name='monitor-xs-ip2',
    version='0.1',
    packages=find_packages(where='src'),
    author='Edoardo Renaldin',
    author_email='edoardo.renaldin@gmail.com',
    package_dir={"": "src"},
    install_requires=[
        'pandas',
        'scipy',
        'numpy',
        'matplotlib',
        'xlsxwriter',
        'periodictable',
        'pysrim'
    ],
    include_package_data=True,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[  # Optional, but a good practice for packaging
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  # Or the license you're using
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',  # Minimum Python version required
)
