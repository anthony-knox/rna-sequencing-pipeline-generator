
Press `esc` and then type in `:wq` and press enter. 
This should exit the vi editor. Type in these commands one by one:

	chmod u+x fastq_dump_Accession.sh
	
	./fastq_dump_Accession.sh @@@ACCESSION_SAMPLES_ABSOLUTEPATH@@@

Your fastq files are now downloading. This may take a long time depending on the size of the files and the number of fastq files you are downloading. While downloading, make sure your computer does not go to sleep as it will interrupt the download process. Setting a screen-saver is OK, but make sure your hard disk does not go to sleep. You may get error messages specifying "timeout exhausted" but as long as there is an output at the end of each sample's fastq-dump that says "Read XXXXXXXXX spots for SRRXXXXXXXX Written XXXXXXXXX spots for SRRXXXXXXXX" then the fastq-dump should have been successful. You can verify by navigating to the kallisto directory with each sample's respective directory, which should now contain the fastq file(s) for that sample.

The fastq files will download to `@@@KALLISTO_ABSOLUTEPATH@@@` in their respective directory.

