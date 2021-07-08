mapt2g.py: a python module that converts transcript coordinates to genomic coordinates, given two 
input files (from problem description):

1. A four column (tab-separated) file containing the transcripts. The first column is the transcript
name, and the remaining three columns indicate it’s genomic mapping: chromosome name, 0-based 
starting position on the chromosome, and CIGAR string indicating the mapping.
2. A two column (tab-separated) file indicating a set of queries. The first column is a transcript
name, and the second column is a 0-based transcript coordinate.

Test input files and output obtained are attached with the solution.

To run the module: download mapt2g.py, input1.txt and input2.txt into working directory first.

```
python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install pandas
python3 mapt2g.py -i1 input1.txt -i2 input2.txt -o output.tsv
```

Assumptions:
* inputs are as indicated, without errors in # columns, or datatype mismatch
* deletions within transcripts are not to be queried
* cigar strings are valid and will only encompass deletions, insertions and (mis)matches
* transcript names are unique

High level description:
* the module is written to store every valid transcript coordinate against a genomic coordinate
* implentation uses python dictionaries, which really are hashtables underneath

Strengths:
* hash-table lookups are pretty quick, so the solution can be engineered towards a one-time build
of the ground truth in input file 1, and subsequent lookups will be against the hash table

Weaknesses:
* this approach takes up space ~ O(entire alignment length). For long transcripts, the hashmap 
building step may (will?) take a while.

Opportunities to improve design/code (high level directions):
* each of the Ms, Is and Ds in the cigar can be implemented like a run-length encoding object, which
can be thought of as a dictionary map for the nearest upstream (5') coordinate + offset. Therefore,
a cigar of 100M starting at position 50 will only occupy one entry in the dictionary (as opposed to
100 entries). Each of the 100 positions can potentially be a binary tree, so after the first lookup
to the nearest genomic coord, the offset will be a log2(offset_length) time lookup.
* ideally there would be classes for transcript entities, and methods written here will be easier to
reason about within that context. 
* unit testing each function reasonably - unfortunately I couldn't do it justice given the time-frame


