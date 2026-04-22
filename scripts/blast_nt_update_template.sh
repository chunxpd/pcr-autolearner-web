#!/bin/sh 
echo "=============================START=============================================="
MOTHDIR=HERE_YOUR_NT_DIRECTORY_PATH
cd $MOTHDIR
cd ori
date
echo "current directory: `pwd`"

# UPBL=`which update_blastdb.pl`
UPBL=/home/i-gene/.local/bin/update_blastdb.pl
$UPBL nt

for tarfile in *.tar.gz; do 
	echo $tarfile
	tar xzf $tarfile
	for suf in nhd nhi nhr nin nnd nni nog nsq; do
		mv ${tarfile%.tar.gz}.${suf} ../
	done
done

mv nt.nal nt.ndb nt.nos nt.not nt.ntf nt.nto taxdb.* ..

cd ../
date
echo "======================END=========================="
exit 0
