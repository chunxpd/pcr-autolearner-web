#!/bin/sh
FLASK_PATH=/home/i-gene/git/ctlabs-primerwrapper
cd $FLASK_PATH
pwd
wget -N -P static/data/ori https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
cd static/data/ori
tar xzf new_taxdump.tar.gz
cd -
PYTHONPATH=. python3 src/taxid_generator.py
