from setuptools import setup, find_packages

#with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

setup(
    name='PTMtools',
    version='0.0.1',
    author='Elly Poretsky',
    author_email='elly.poretsky@gmail.com',
    #packages=['PTMtools'],
    description='Toolking for handling PTM data',
    #long_description="long_description",
    #long_description_content_type="text/markdown",
    url='https://github.com/eporetsky/ptmtools',
    project_urls = {
        "Bug Tracker": "https://github.com/eporetsky/ptmtools/issues"
    },
    license='GNU3',
    
    install_requires=[
        'pandas',
        'numpy',
        'biopython',
        'h5py',
        'seaborn',
        'matplotlib',
    ],


    entry_points={
        'console_scripts': ['ptmtools=functions:load_fasta'],
        #'console_scripts': ['ptmtools load_fasta=functions:load_fasta']
    }
)