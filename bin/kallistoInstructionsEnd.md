Press `esc` and then type in `:wq` and press enter. 
This should exit the vi editor. Type in these commands one by one:

	chmod u+x kallisto.sh

	./kallisto.sh @@@KALLISTO_SAMPLES_ABSOLUTEPATH@@@

The aligned reads are now placed in a subdirectory named `output` within each sample's respective directory in `@@@KALLISTO_ABSOLUTEPATH@@@`

If you would like to see a combined text file of the aligned reads, enter the following commands into terminal/command prompt:

	mkdir @@@KALLISTO_ABSOLUTEPATH@@@/Abundance_files_named
	
	cd @@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@
	
	vi rename_move_abundance_files.sh

This brings you to the vi editor to create a shell script. Press `i` then copy and paste the code below:

	#!/bin/bash
	while IFS= read -r line || [[ -n "$line" ]]; do
		cd @@@KALLISTO_ABSOLUTEPATH@@@/$line/output
		cp abundance.tsv @@@KALLISTO_ABSOLUTEPATH@@@/Abundance_files_named/$line.tsv
	done < "$1"

Press `esc` and then type in `:wq` and press enter. 
This should exit the vi editor. Type in these commands one by one:

	chmod u+x rename_move_abundance_files.sh
	
	./rename_move_abundance_files.sh @@@KALLISTO_SAMPLES_ABSOLUTEPATH@@@
	
	cd @@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@
	
	vi merge_tsv_files.sh

This brings you to the vi editor to create a shell script. Press `i` then copy and paste the code below:

	#!/bin/bash
	OUT=@@@KALLISTO_ABSOLUTEPATH@@@/Abundance_files_named/all_merged.csv
	touch $OUT
	
	for file in @@@KALLISTO_ABSOLUTEPATH@@@/Abundance_files_named/*.tsv
	do
	        filename=$(basename "$file")
	        fname="${filename%.*}"
	        echo $fname
	        paste $OUT <(awk -v OFS=',' '{print $1,$2,$3,$4,$5,""}' $file) >> $OUT.tmp
	        mv $OUT.tmp $OUT
	done

Press `esc` and then type in `:wq` and press enter. 
This should exit the vi editor. Type in these commands one by one:

	./merge_tsv_files.sh

You now have a csv file that you can use to analyze the raw data between samples via Excel. Each sample's data is placed horizontally to the next sample. Note that the data is placed horizontally based on lexicographical order of your samples (that which shows up in your terminal output after running the script). If you only want specific columns, you can edit the shell script at the line that begins with "paste ..." and only include the column numbers that you want. For example, if you only want the tpm column, include only $5. Be careful in drawing conclusions about the data because these are raw counts that have not been normalized or applied to statistical models that are used in the DESeq2 package.

