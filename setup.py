import distutils
from distutils.core import setup

# The main call
setup(name='qatoolkit',
      version ='0.1.2',
      license = "GPL",
      description = "A set of handy Python utility scripts for DESDM",
      author = "Felipe Menanteau",
      author_email = "felipe@illinois.edu",
      packages = ['qatoolkit'],
      package_dir = {'': 'python'},
      scripts = ['bin/color_tile',
                 'bin/compare_corners_NewFramework',
                 'bin/display_DECam_MEFs',
                 'bin/display_DECam_detection_NOMAD',
                 'bin/make_png_focus_chips_SQL',
                 'bin/make_png_focus_chips_bypath',
                 'bin/projectDECamPNG',
                 'bin/projectDECamPNG_SQL',
                 'bin/projectDECamPNG_bypath'],
      data_files=[('ups',['ups/qatoolkit.table']),
                  ('etc',['etc/default.stiff','etc/default.swarp'])]
      )           
                 

