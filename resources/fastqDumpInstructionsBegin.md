#Fastq Dump Instructions
Please visit [the SRA toolkit download page](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) to download the latest version of the SRA toolkit.
Unzip the file, rename the main directory to `sratoolkit` and move it to `@@@PWD_ABSOLUTEPATH@@@`

Make sure that the `samples_Accession.txt` file within the `kallisto` directory contains the accession numbers for each sample. 
If it does not, simply re-enter your sample names with correct formatting and click the `Enter Sample Names` button again. Note that there may be incorrectly named directories created previously if you change your sample names, so you will have to remove those yourself.

Open a new terminal/command prompt window and enter these commands, one by one:

	cd @@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@
	
	PATH=$PATH:@@@PWD_ABSOLUTEPATH@@@/sratoolkit/bin
	
	vi fastq_dump_Accession.sh

This brings you to the vi editor to create a shell script. Press `i` then copy and paste the code below:

