SAP 
=====

A DNA sequence analysing & mutation predicting tool.


How to use
=====

Compile
-----

    make

Basic usage
-----

1. Transform Fasta format into FDA format.
    

    FastaToFDA REFERENCE_INPUT_FILE_NAME > REFERENCE_FDA.fda


2. Transform Fastq format into FDQ format.

    
    FastqToFDQ INPUT_FILE_NAME > INPUT_FDQ.fdq


3. Run SAP Mapper.

    
    Mapper -i INPUT_FDQ.fdq -r REFERENCE_FDA.fda -o RESULT.txt


