# RNA-Sequencing Pipeline Generator

The RNA-Sequencing Pipeline Generator is a program that provides detailed instructions for setting up an analytical platform for RNA-Sequencing analysis, utilizing the Kallisto and DESeq2 packages. This program is designed for researchers who are relatively new to RNA-Sequencing analysis and who have little to no experience with the Linux/Unix or R programming languages. Code will be provided that is specific to your own file structure, allowing for a seamless and streamlined process that includes downloading fastq files from NCBI using the SRA toolkit (if necessary), performing sequence alignment with the Kallisto package, then performing differential expression analysis with the DESeq2 package. 

Detailed instructions are provided to the user for each required step. Additionally, shell scripts are provided to automate various steps, such as downloading fastq files for a list of specified NCBI accession numbers. Output from this program includes transcript level abundances, a normalized counts matrix, a list of up- and down-regulated genes, hierarchical clustering heatmaps, PCA plots, and bi-modality analysis via the SIBER and DEXUS packages. Sample code is provided for each of these types of analyses, but the user can explore each function further to specify parameters that meet their own needs. Additionally, further resources are provided which allow for downstream analysis on your normalized counts matrix.

## Getting Started

1. Download the JAR file containing the executable program [here](https://github.com/anthony-knox/rna-sequencing-pipeline-generator/blob/master/RNA-Sequencing_Pipeline_Generator.jar).
2. Either double click the file or `cd` to the directory that contains the JAR file and enter `java -jar RNA-Sequencing_Pipeline_Generator.jar`

### Prerequisites

#### System Requirements

1. Installation of Java 8 or later JRE
	* [Java 8 JRE](https://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html)
	* To check what version of Java you have, open a terminal or command prompt window and enter ```java -version```
2. Enough disk space to store fastq files (usually ~1 GB per sample)
3. Markdown (optional)
	* The output pipeline file is easier to read by converting the text file to HTML
	* Visit [What Is an MD File?](https://www.lifewire.com/md-file-4143558) to learn more
	* If using Google Chrome, the [Markdown Preview Extension](https://chrome.google.com/webstore/detail/markdown-preview/jmchmkecamhbiokiopfpnfgbidieafmd?hl=en) is very useful to view markdown files
		* Download the extension, then navigate to `chrome://extensions` in your browser, click `Details`, and check `Allow access to file URLs`

_Currently this program is only supported on Mac OS. Support for Windows will be provided in the future. Windows users will need to look into installing [Cygwin](https://www.howtogeek.com/howto/41382/how-to-use-linux-commands-in-windows-with-cygwin/) and will need to alter the bash scripts to support Windows functionality._

#### Files and Directories Requirements

1. The absolute path of the present working directory
	* Example: `/Users/userName/Desktop/RNA_Seq_Analysis`
	* If using either Windows or Mac, enter the path with forward slashes as shown
	* Initially creating an empty directory is suggested
2. The sample accession numbers (if downloading from NCBI), sample names, and conditions for each sample, separated by a comma: `<accession#>,<samplename>,<condition>`
	* Example: `SRR1027605,C9_52i_run1,exp`
3. If downloading from NCBI is not necessary, the fastq files must be gzipped and named according to the sample names
	* Example: `C9_52i_run1.fastq.gz`
4. Whether your reads are single-end or paired-end

All other necessary components will be downloaded by following the instructions provided in the program.

### Usage

As an example of running this program, you can use publicly available single-end RNA-Sequencing data from iPSC derived motor neuron cultures from C9ORF72 carriers. More information about this project can be found [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA227155). 

To begin, run the JAR file and enter the absolute path of your present working directory (see Prerequisites: Files and Directories #1 above). Then click `Create File Structure`

Provided below are the sample accession numbers, sample names, and conditions, which is used as input in the next section. Copy and paste these lines into the text area and click `Enter Sample Names`.

SRR1027605,C9\_52i\_run1,exp  
SRR1027606,C9\_52i\_run2,exp  
SRR1027603,C9\_30i\_run1,exp  
SRR1027604,C9\_30i\_run2,exp  
SRR1027601,C9\_29i\_run1,exp  
SRR1027602,C9\_29i\_run2,exp  
SRR1027599,C9\_28i\_run1,exp  
SRR1027600,C9\_28i\_run2,exp  
SRR1027597,Control\_83i\_run1,ctl  
SRR1027598,Control\_83i\_run2,ctl  
SRR1027595,Control\_14i\_run1,ctl  
SRR1027595,Control\_14i\_run2,ctl  
SRR1027593,Control\_03i\_run1,ctl  
SRR1027594,Control\_03i\_run2,ctl  
SRR1027591,Control\_00i\_run1,ctl  
SRR1027592,Control\_00i\_run2,ctl  

In the following section check `I need to download my .fastq files from NCBI` and `Single-End` and click `Generate fast dump instructions`. Note that fastq files are relatively large files and the fastq-dump will take a decent amount of time. You can leave this fastq-dump script running overnight and it should be finished by the following morning.

Further detailed instructions are provided while progressing through the program. If you require any assistance, please don't hesitate to contact the repository owner by email. Providing a detailed question will help with any troubleshooting.

## Contributing

Pull requests are welcome. Please first discuss the change you wish to make via email with the owner of this repository before making a change.

Please follow the code of conduct available from the [Contributor Covenant](https://www.contributor-covenant.org/version/1/4/code-of-conduct.html)

## Authors

* **Anthony Knox** - *Initial work* - [https://github.com/anthony-knox](https://github.com/anthony-knox)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

