## Navigation ##
  * [Return to Table of Contents](TableOfContents.md)
  * [Jump to Understanding Kojak Output](KojakOutput.md)

## Introduction ##

Kojak uses the Percolator algorithm to validate peptide spectrum matches and assign a q-value to all search identifications. The q-value can be used to apply a false discovery rate (FDR) cutoff. Percolator is not maintained on this site. It can be obtained for free at: http://sourceforge.net/projects/percolator/files/

Additional Percolator information is found at: http://per-colator.com/

**SPECIAL NOTE: Kojak now supports formatting Percolator input to work with different versions of Percolator. Please set the percolator\_version parameter in your Kojak configuration file to the version of Percolator that you use.**


## Installation ##

Percolator comes as pre-compiled binaries or source code. We recommend using the binaries. Instructions for compiling the souce code are beyond the scope of this document.

## Use ##

Percolator operates from the command line. To facilitate using Percolator, Kojak produces an optional [Percolator-ready output file](KojakOutput.md). The results from Percolator are displayed on the screen and can be captured by piping to a text file.

Detailed usage of Percolator is documented on the [Percolator web page](http://per-colator.com/). Below is the typical usage statements for validating Kojak output.

If you use Percolator version 2.04 and older:

```
percolator -j your_Kojak_percolatorInput_file > your_desired_output_filename.txt
```

For Percolator versions 2.05 and newer:

```
percolator your_Kojak_percolatorInput_file > your_desired_output_filename.txt
```