# Fast Genome Profiler (Fast-GeP)
Fast and easy tool to infer high-resolution genealogical relationship of bacterial isolates from WGS data.

## Introduction

Fast-Genome Profiler (Fast-GeP) is the first allele calling program that using genome-by-genome algorithm to transform bacterial WGS assembly data into allele profiles. 
It is written in Perl and has been tested in MacOS (OS X 10.9.5 and High Sierra v 10.13.1) and Ubuntu (14.04 and 16.04).

## Motivation
Large-scale WGS dataset based studies are becoming increasingly common. WGS based pathogen surveillance and outbreak investigations create a demand for highly discriminative and time-efficient bioinformatics tools to transform large amounts of sequencing data into usable biological information such as relationship of the isolates.
Installation

## Prerequisites
Before start, you need to make sure the following three programs were full functional in your system:
   * [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
   * [DIAMOND](https://github.com/bbuchfink/diamond)
   * [MAFFT](https://mafft.cbrc.jp/alignment/software/)

## Usage
Once the three dependent external tools were installed properly, you could copy the program to a directory in which your the data files stored. 
Initiate the analysis by using the command: 

    perl fast-GeP.pl <options>
    
## Examples:
    perl fast-GeP.pl -g list.fas.txt -r reference.gbk

This commond will build an _ad hoc_ scheme from `reference.gbk` to transform the assembly files in `list.fas.txt` into allele profiles using BLAST+.

## Citation



