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

STEP 1: Transform Fasta format into FDA format.

    FastaToFDA REFERENCE_INPUT_FILE_NAME > REFERENCE_FDA.fda

STEP 2: Transform Fastq format into FDQ format.

    FastqToFDQ INPUT_FILE_NAME > INPUT_FDQ.fdq

STEP 3: Run SAP Mapper.

    Mapper -i INPUT_FDQ.fdq -r REFERENCE_FDA.fda -o RESULT.txt

STEP 4: Run Predictor.
    
    Predictor -i RESULT.txt -o VARIATION.txt -r REFERENCE_FDA.fda
    
The file VARIATION.txt contains the final result.

Options for SAP Mapper
-----

Options for SAP Predictor
-----

Options for SAP SNPFilter
-----
