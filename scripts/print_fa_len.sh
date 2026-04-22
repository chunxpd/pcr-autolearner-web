#!/bin/sh

for in_fa in $*; do
	echo $in_fa
	grep -P "^>" $in_fa | cut -f 1,1 -d " " | sed -e "s/complement(//" -e "s/)//g" | sed -e "s/:/ /" -e "s/-/ /" | gawk '{print $1":"$2"-"$3" "$3-$2}' | sort -k2,2 -nr | head
	#grep -P "^>" $in_fa | cut -f 1,1 -d " " | cut -f 2,2 -d ":" | sed -e "s/complement(//" -e "s/)//g" | gawk -F "-" '{print($2-$1)}' | sort -n -r | head
done
