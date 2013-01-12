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

Suppose the reference file name is REFERENCE_INPUT_FILE_NAME (for example, human gene), 
and INPUT_FILE_NAME contains the reads that needed to be mapped onto the reference.

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

*   -f 

    SAP cuts read into small pieces.
    Usually, one SNP between small piece and reference can be tolerated.
    When FastMap (-f) is enabled, any SNP between small piece and reference will not be tolerated, 
    which can greatly accelerate the mapping process, and reduce the coverage of mapping.

*   -t THREAD_COUNT

    The number of threads when mapping.
   
*   -H HASH_SIZE
    
    The size of hash table.
    Can be any number between 20 and 30.
    The actual size of hash table will be 2^HASH_SIZE.
    Larger size of hash table will lead to faster mapping when reference is large.

*   -C CUT_COUNT

    The number of pieces that every read will be cut into.
    Each piece will be looked up in the hash table.
    Usually, larger number of pieces will lead to higher coverage of mapping.

Options for SAP Predictor
-----

Options for SAP SNPFilter
-----
