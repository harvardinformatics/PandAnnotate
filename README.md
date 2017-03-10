# PandAnnotate
A python-based alternative to Trinotate for merging herterogeneous transcriptome assembly annotation tables

##Background
Integrating information across annotations and nucleotide/protein database searches is a key step to understanding the functional content of a de novo transcriptome assembly, identifying contigs that may come from exogenous sources (e.g. bacteria, viruses, parasites, symbionts), establishing criteria to remove contigs from an assembly, and for performing many downstream analyses. The developers of Trinity created Trinotate, which adds annotations to a SQLite database, that is then written to, wait,....an Excel file. While Trinotate is relatively easy to use, it compresses many fields of a database search into a single column, delimited by heterogeneous separators that are difficult to parse, let alone load into R. Furthermore, it lacks flexibility by requiring the use of sequence-databases customized for Trinotate, not allowing integration of blastn searches, or other arbitrarily structured annotation information such as expression data. The goal of PandAnnotate is to provide a more flexible, pythonic implementation of database/search merges, and to write the merged database to a text file that is easy to load into R for further querying and analysis.

##Command line arguments for PandAnnotate.py are as follows: 

* **-f,--transcriptome_fasta**	fasta file of assembled contigs

* **-c,--control_file**		control file (see explanation below)

* **-o,--outtable**		name to which merged annotation table will be written

# Building a control file
The control file consists of four tab-separated columns per annotation table that is provided to PandAnnotate.py:

* full path to file

* search type: current options are blastn, blastx, blastp, transdecoder, and custom

* prefix: a string you want pre-pended to the column names (either the defaults for blast,pfam and transdecoder, or for a supplied set of column names). Since pfam and transdecoder are only run once--such that there will not be redundant column names--one can omit a prefix for these searches. However, a prefix is necessary if you are integrating multiple blast searches, in order to prevent redundant column names and you confusing one blast search with another.

* column names: for blast,pfam, and transdecoder, this can be left blank,as there are default column labels. Otherwise, this is a comma-separated list of column names in the order they appear in the table.

**NOTE:** blast search tables need to be in outfmt6 format. Column labels must be provided if you use a custom-defined set of reported fields rather than the outfmt6 defaults. 		

## Coming soon:

* Automatic pulldown of GO terms from remote databases
* Functions for annotation table filtering
* Output of table that lists contigs and their pass/fail status with respect to filtering
* Descriptive statistics
