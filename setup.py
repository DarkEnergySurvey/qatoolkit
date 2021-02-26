import distutils
from distutils.core import setup
import glob

bin_files = glob.glob("bin/*")
# The main call
setup(name='qatoolkit',
      version ='3.0.0',
      license = "GPL",
      description = "A set of handy Python utility scripts for DESDM",
      author = "Robert Gruendl",
      author_email = "gruendl@illinois.edu",
      packages = ['qatoolkit'],
      package_dir = {'': 'python'},
      scripts = bin_files,
      data_files=[('ups',['ups/qatoolkit.table'])]
      )
