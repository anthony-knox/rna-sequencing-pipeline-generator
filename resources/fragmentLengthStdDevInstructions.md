#Instructions for fragment length and standard deviation

Since your fastq files are single end, you will need to specify your estimated average fragment length and estimated standard deviation of fragment length.
If you don't know these values, open a terminal/command prompt window and `cd` to one of your sample directories in `@@@KALLISTO_ABSOLUTEPATH@@@`.
Unzip all sample files so that they have a .fastq extension. Then copy and paste this command in terminal/command prompt (this will take a few minutes complete):

	for F in *.fastq
	do 
	echo  "$F   "
	awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sq/n-m*m);}' $F
	done

You can repeat this for multiple sample files if you want to get more accurate values for your average fragment length and standard deviation.
Your fragment length average should be somewhere around 50 - 200, and your standard deviation should be somewhere around 0 - 20.

