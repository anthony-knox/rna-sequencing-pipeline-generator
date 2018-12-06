#Instructions to move fastq files to proper directory

Make sure that your sample names in your `samples_kallisto.txt` file match your file names.
In the text file, your sample names will only include the basename without file extensions (without .fastq.gz). For example, a file name should be `sampleName.fastq.gz` for single-end data or `sampleName_1.fastq.gz` and `sampleName_2.fastq.gz` for paired-end data (gzip your files if they are not already).

Place all fastq files in: 
`@@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@`

Open a new terminal/command prompt window and enter these commands, one by one:

	cd @@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@
	
	vi mv_fastq_to_sample_dir.sh

This brings you to the vi editor to create a shell script. Press `i` then copy and paste the code below:

	#!/bin/bash
	while IFS= read -r line || [[ -n "$line" ]]; do
	        echo $line
	        mv @@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@/${line}* @@@KALLISTO_ABSOLUTEPATH@@@/$line
	done < "$1"

Press `esc` and then type in `:wq` and press enter. 
This should exit the vi editor. Type in these commands one by one:

	chmod u+x mv_fastq_to_sample_dir.sh
	
	./mv_fastq_to_sample_dir.sh @@@KALLISTO_SAMPLES_ABSOLUTEPATH@@@

Your fastq files should now be placed in their respective directories within the kallisto directory.

