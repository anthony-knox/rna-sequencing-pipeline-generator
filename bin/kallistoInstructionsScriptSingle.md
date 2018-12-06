	#!/bin/bash
	while IFS= read -r line || [[ -n "$line" ]]; do
		echo $line
		cd @@@KALLISTO_ABSOLUTEPATH@@@/$line
		kallisto quant -i @@@KALLISTO_ABSOLUTEPATH@@@/ensemblIndex.idx -o output --single -l @@@FRAGMENT_LENGTH@@@ -s @@@STD_DEV@@@ ${line}.fastq.gz
	done < "$1"

