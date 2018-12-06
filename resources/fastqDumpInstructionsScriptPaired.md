	#!/bin/bash
	ulimit -n 10000
	while IFS= read -r line || [[ -n "$line" ]]; do
	        ACCESSION=$(echo "$line" | cut -d',' -f1)
	        NAME=$(echo "$line" | cut -d',' -f2)
	
	        echo "Downloading $ACCESSION for sample $NAME"
	        fastq-dump -I --split-files --gzip -O @@@KALLISTO_ABSOLUTEPATH@@@/$NAME $ACCESSION
	        mv @@@KALLISTO_ABSOLUTEPATH@@@/$NAME/${ACCESSION}_1.fastq.gz @@@KALLISTO_ABSOLUTEPATH@@@/$NAME/${NAME}_1.fastq.gz
	        mv @@@KALLISTO_ABSOLUTEPATH@@@/$NAME/${ACCESSION}_2.fastq.gz @@@KALLISTO_ABSOLUTEPATH@@@/$NAME/${NAME}_2.fastq.gz
	done < "$1"

