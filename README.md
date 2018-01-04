# Fast Genome Profiler (Fast-GeP)
Fast and easy way to infer high-resolution genealogical relationship of bacterial isolates with WGS assembly data.

## Introduction

Fast-Genome Profiler (Fast-GeP) is the first genome-by-genome allele-calling program that can transform bacterial WGS assembly data into allele profiles. 
It is written in Perl and has been tested in MacOS (OS X 10.9.5 and High Sierra v 10.13.1) and Ubuntu (14.04 and 16.04).

## Motivation
Large-scale WGS dataset based studies are becoming increasingly common. WGS based pathogen surveillance and outbreak investigations create a demand for highly discriminative and time-efficient bioinformatics tools to transform large amounts of sequencing data into usable biological information such as relationship of the isolates.

## Prerequisites
Before start, you need to make sure the following three programs were full functional in your system:
   * [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
   * [DIAMOND](https://github.com/bbuchfink/diamond)
   * [MAFFT](https://mafft.cbrc.jp/alignment/software/)

## Usage
Let's assume you have put the `fast-GeP.pl` file in your PATH. If not, or if you prefer, you could always put the `fast-GeP.pl` file along with your other input files, and use commands like:

    perl fast-GeP.pl -g list.fas.txt -r reference.gbk

Here are some example commands:

Run the analysis using BLAST+ as aligner (_i.e._ build an _ad hoc_ scheme from `reference.gbk` to transform the assembly files in `list.fas.txt` into allele profiles using BLAST+):

    fast-GeP -g list.fas.txt -r reference.gbk

Run the analysis using DIAMOND as aligner:

    fast-GeP -g list.fas.txt -r reference.gbk -d

Run the analysis using BLAST+ as aligner and produce a 'long' results:

    fast-GeP -g list.fas.txt -r reference.gbk -l

Run the analysis using BLAST+ as aligner and do not produce pariewise comparison files:

    fast-GeP -g list.fas.txt -r reference.gbk -n

## Citation
_Genome-by-genome approach for fast bacterial genealogical relationship evaluation_
