	#!/bin/bash
	ulimit -n 10000
	while IFS= read -r line || [[ -n "$line" ]]; do
		ACCESSION=$(echo "$line" | cut -d',' -f1)
		NAME=$(echo "$line" | cut -d',' -f2)
	
		echo "Downloading $ACCESSION for sample $NAME"
		fastq-dump -I --gzip -O @@@KALLISTO_ABSOLUTEPATH@@@/$NAME $ACCESSION
		mv @@@KALLISTO_ABSOLUTEPATH@@@/$NAME/$ACCESSION.fastq.gz @@@KALLISTO_ABSOLUTEPATH@@@/$NAME/$NAME.fastq.gz
	done < "$1"

