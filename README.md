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
    FastMap mode.  
    SAP cuts read into small pieces.
    Usually, one SNP between small piece and reference can be tolerated.
    When FastMap (-f) is enabled, any SNP between small piece and reference will not be tolerated,
    which can greatly accelerate the mapping process, and reduce the coverage of mapping.


*   -t THREAD_COUNT  
    The number of threads when mapping.


*   -H HASH_SIZE  
    The size of hash table, which can be any number between 20 and 30.  
    The actual size of hash table will be 2^HASH_SIZE.
    Larger size of hash table will lead to faster mapping when reference is large.


*   -C CUT_COUNT  
    The number of pieces that every read will be cut into.  
    Each piece will be looked up in the hash table.
    Usually, larger number of pieces will lead to higher coverage of mapping.


*   -G GAP_RATIO  
    Maximum gap ratio.  
    The maximum percentage of gaps in a single read that can be tolerated by SAP Mapper.
    The gap ratio is defined as MaxGapLength/ReadLength.


*   -i FILE_NAME  
    Input file name (the file which contains reads, should be in FDQ format).


*   -r FILE_NAME  
    Reference file name (the file which contains reference, should be in FDA format).


*   -o FILE_NAME  
    Output file name.


*   -p PIECE_SIZE  
    The size of small pieces when mapping.  
    Smaller size of pieces will lead to slower mapping and higher coverage.


*   -h  
    Help.


Options for SAP Predictor
-----

*   -i FILE_NAME  
    Input file name (the output of SAP Mapper).


*   -r FILE_NAME  
    Reference file name (the file which contains reference, should be in FDA format).


*   -o FILE_NAME  
    Output file name.


*   -Q QUALITY  
    Minimum quality to validate a read.


*   -m READ_DEPTH  
    Mimimum depth of reads to validate a insertion/deletion/SNP.  
    Here the depth of reads means the number of reads that covers the location of insertion/deletion/SNP.


*   -h  
    Help.

Options for SAP SNPFilter
-----

*   -i FILE_NAME  
    Input file name (the output of SAP Predictor).


*   -o FILE_NAME  
    Output file name.


*   -m READ_DEPTH  
    Minimum depth of reads to call a SNP.


*   -M READ_DEPTH  
    Maximum depth of reads to call a SNP.


*   -s SCORE  
    Maximum score to call a SNP.  
    For the output of SAP Predictor, higher score means lower possiblity that a SAP happens on the position.


*   -h  
    Help.
