# fusia
Find Unique Substrings In Assemblies

## Authors
Sarah Spencer <sjspen@gmail.com>, Sept. 2019

## Quickstart
Check for uniqueness of an input sequences within a collection of fasta files:
~~~
$ fusia.py substring -i <INDIR> -s <SEQUENCE> -o <OUTFILE>
~~~

Find unique substrings across a collection of fasta files:
~~~
$ fusia.py unique -i <INDIR> -o <OUTDIR>
~~~

## Dependencies
* Python v3.7.4
* Pandas v0.25.1
* Biopython v1.74
* Mauve v2.4.0

