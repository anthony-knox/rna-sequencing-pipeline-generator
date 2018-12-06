#Instructions for Kallisto alignment

[Kallisto manual](https://pachterlab.github.io/kallisto/manual)

Please visit [the Kallisto download page](https://pachterlab.github.io/kallisto/download) to download Kallisto on your specific platform.

Visit [the Ensembl FTP download page](http://uswest.ensembl.org/info/data/ftp/index.html) and choose a reference genome to download from the cDNA column (the FASTA link). You will download the one that contains all genes (For example: Homo_sapiens.GRCh38.cdna.all.fa.gz).
Once the file is downloaded, move it to `@@@KALLISTO_ABSOLUTEPATH@@@` and rename it to `ensemblRefGenome.fa.gz`. If you obtain your reference genome elsewhere, you can name it whatever you would like. Just be sure to change the later code to match your name.

Open a new terminal/command prompt window and enter these commands, one by one:

	cd @@@KALLISTO_ABSOLUTEPATH@@@
	
	kallisto index -i ensemblIndex.idx ensemblRefGenome.fa.gz

After your index file is created, you are ready to perform the quantification algorithm. Enter these commands, one by one:

	cd @@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@

	vi kallisto.sh

This brings you to the vi editor to create a shell script. Press `i` then copy and paste the code below:

