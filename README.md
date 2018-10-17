# Demultiplexing

An index tag, also known as a barcode, is a unique string of 6-12 nucleotides added to each sample sequence to allow different samples in lanes of an Ilumina flow cell to be identified. This facilitates multiplexing, in which large numbers of libraries are pooled
and sequenced simultaneously during a single run. Demultiplexing is the reverse process of multiplexing, in which
sequence reads are divided into separate files based on their index tag. This process is imperfect; often times moderate amounts of sequencing data must be discarded because they fail to meet quality checks. 

This program is designed to accept paired end reads and check if they have been properly mapped. It requires two gzipped fastq files of biological reads and two gzipped fastq files of index tag reads. The program also requires an index library with the indexes stored in the 5th column and a qscore of the user's preference. If they indexes are properly mapped to their expected sequence, the index ID will be recorded in the header line of the biological read and the forward and reverse biological fastq files will get written to a gzipped file of their respective index. 

This program is filters out reads for the following reasons:
  - There is an undertermined (N) nucleotide in one of the index sequences
  - The mean quality score one of the indexes is below the user-designated quality score threshold
  - At least one of the indexes does not appear in the index library
  - Index hoppping between one or both indexes

Whenever a read is filtered out for one of the above phenomena, that reason will be commented on the third line of its fastq file. More than one reason can be recorded per file, as a line must pass through each quality filter. After each reason for failing a quality check has been recorded, the entire forward and reverse file will be written to their respective "bad" gzipped file. As with the good files, the index sequence will be appended to the header line of the biological read.
  
The program will return the following reports:
  - The raw count of each index as well as its relative percentage of the total amount of lines seen
  - The proportion of the bad files expressed as a percentage
  - The number of instances a file was filtered for having an N, a qscore below the threshold, illigitimate indexes, or index hopping
  
Note that only a single instance is recorded per filtering event. Also of note, the counts that are returned for reads that failed tests are not mutually exclusive (e.g. a file could with a low quality score could have also experienced index hopping). To quantify how many times a specific quality failure or combinations of quality failures occured, you can use bash commands to extract the phrases "Undetermined", "Low Quality", "Index not recognized", or "Index hopping" from the third line of either "bad" file.
