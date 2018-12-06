	#!/bin/bash
	while IFS= read -r line || [[ -n "$line" ]]; do
		echo $line
		cd @@@KALLISTO_ABSOLUTEPATH@@@/$line
		kallisto quant -i @@@KALLISTO_ABSOLUTEPATH@@@/ensemblIndex.idx -o output ${line}_1.fastq.gz ${line}_2.fastq.gz
	done < "$1"

\*Note that in the script above, your paired-end fastq files should be named as `sampleName_1.fastq.gz` and `sampleName_2.fastq.gz`

