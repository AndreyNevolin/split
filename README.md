#File splitting utility


## Description
Utility that splits a file of "elements" of any type into pieces of roughly
equivalent sizes. Each piece would contain an integer number of "elements"
(i.e. one single element couldn't be divided between pieces. Each element is
required to belong to some piece as a whole).

Exact definition of "elements" (and thus complete purpose of the tool) is
provided by a party that builds the tool.

The tool reads data from input file and writes data to output files in aligned
chunks of fixed size (except maybe the last chunk read from/appended to a file).
The size of a chunk can be specified by a user. By default it's 4Mb.

The tool works best when the chunk size is much bigger than the size of any
"element".

## Configuring
The tool is incomplete. It's a prototype which is agnostic of exact source file
format (i.e. agnostic of "elements" definition). To make the tool complete, one
needs to implement "split_FindBound()" function which is defined in
"find_bound.cpp" file. This file comes with a reference implementation of
"split_FindBound()" intended for splitting files in FASTA format (bioinformatics
format for representing nucleotide or peptide sequences). To build the tool for
splitting files of different format, one needs to replace the contents of
"split_FindBound()" with corresponding code.

The tool reads data from input file and writes data to output files in aligned
chunks of fixed size (except maybe the last chunk read from/appended to a file).
The size of a chunk can be specified by a user. By default it's 4Mb.

To ensure that output files contain an integer number of "elements" each, the
tool somehow needs to recognize bounds of individual elements. When last chunk of
data is added to an output file, projected bound of a file might be shifted up or
down to make the file include integer number of elements. That's where
"split_FindBound()" comes to play. Given a buffer with data and projected output
file bound, it should be able to find an element bound which is close to the
projected bound.

For more information about "split_FindBound()" see comments inside "find_bound.cpp

## Building
There are two options:
1) run ```make``` to build release version of the tool
2) run ```make debug``` to build debug version of the tool

The debug version comes with a symbol table and lots of internal sanity checks

## Using the tool
```
split -n <number of pieces> [-od <output directory>] [-of <basis for output file name>] [-cs <chunk size>] <path to file to split>

OPTIONS:
   -n          Number of pieces to produce. Each piece will be placed into
               a separate file named "<file name>.<number>", where "<file
               name>" is a name provided through "--of" option or the
               name of the input file (if "--of" option is not used).
               "<number>" is a sequential number of a piece
       --od    Path to output directory. By default current directory will
               be used for output
       --of    Basis for output file names. Output files will be named
               "<file name>.<number>", where "<file name>" is a string
               provided through this option or the name of the input file
               (if this option is not used)
       --cs    Data chunk size. Data will be read from input file and written
               to output files by chunks of this size (except maybe the last
               chunk read from/written to file). The size may be provided in
               different units: 512B, 4K, 8M, 1G ("B" for bytes, "K" for
               kilobytes, "M" for megabytes, and "G" for gigabytes). If
               units identifier is omitted, byte units are implied. The
               default value for this option is 4M
```

## Copyright
Copyright © 2016 Andrey Nevolin, https://github.com/AndreyNevolin
 * Twitter: @Andrey_Nevolin
 * LinkedIn: https://www.linkedin.com/in/andrey-nevolin-76387328