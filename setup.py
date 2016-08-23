import distutils
from distutils.core import setup

# The main call
setup(name='qatoolkit',
      version ='0.2.6',
      license = "GPL",
      description = "A set of handy Python utility scripts for DESDM",
      author = "Robert Gruendl",
      author_email = "gruendl@illinois.edu",
      packages = ['qatoolkit'],
      package_dir = {'': 'python'},
      scripts = ['bin/assess_SE_products.py',
                 'bin/check_flats.py',
                 'bin/flats_analysis.py'],
      data_files=[('ups',['ups/qatoolkit.table'])]          
      )     
