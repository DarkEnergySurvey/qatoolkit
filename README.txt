To set up for testing:
----------------------
 setup -r ~/DESDM-Code/devel/qatoolkit/trunk
   or
 setup -r ~/DESDM-Code/devel/qatoolkit/tags/XX.YY.ZZ

To install:
-----------
 
 python setup install --home=$HOME/Python
  or
 python setup.py install --prefix=$PRODUCT_DIR --install-lib=$PRODUCT_DIR/python 


Scripts: (inside bin/)
----------------------

 - assess_SE_products.py : Robert's assesments script. Perfoms an
 assessent of exposures from a first/final cut run.  The quality of
 the exposures is based upon the seeing (FWHM), background, and
 extinction due to clouds


