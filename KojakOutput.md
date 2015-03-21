# Interpreting Kojak Results #

## Navigation ##
  * [Return to Table of Contents](TableOfContents.md)


## Overview ##

Kojak produces two types of output files: 1) a summarized list of all spectra results, and 2) a set of PSM summary files formatted for further analysis using [Percolator](Percolator.md). Both output file types are presented in ASCII text and are tab-delimited. They can alternately be opened in most spreadsheet software (such as Microsoft Excel) for easy viewing.

## Output File Type #1: Kojak PSMs and Scores ##

The Kojak summary file (typically appended with the text "-kojak") contains the results for all scans searched. If the fields after the scan number contain zeroes, this means that the spectrum had no PSMs of any sort. This could happen if the database does not contain any peptides of the expected mass. For the remainder of the spectra, the results are as follows:

  1. **Scan Number** - The scan even number.
  1. **Obs Mass** - The observed precursor mass for the spectrum.
  1. **Charge** - The observed precursor charge for the spectrum.
  1. **PSM Mass** - The calculated exact mass of the matching peptide(s).
  1. **PPM Error** - The mass error between Obs Mass and PSM Mass, in parts per million.
  1. **Score** - The Kojak cross correlation score. For cross-linked PSMs, this is the summed Score of each peptide.
  1. **dScore** - The delta score, or difference between the Score of this PSM and the Score of the next best PSM.
  1. **Pep. Diff.** - For cross-linked PSMs, the difference between the Score and the contribution from the best scoring peptide of a cross-link, i.e. the Score of the lower scoring peptide.
  1. **Peptide 1** - The longer peptide sequence in a cross-link, or the only peptide sequence of all other PSMs.
  1. **Link 1** - The site of the cross-link, relative to the peptide sequence, or -1 if there is no cross-link.
  1. **Protein 1** - The protein(s) that contain Peptide 1.
  1. **Peptide 2** - The shorter peptide sequence in a cross-link.
  1. **Link 2** - The site of the cross-link, relative to the peptide sequence.
  1. **Protein 2** - The protein(s) that contain Peptide 2.
  1. **Linker Mass** - The mass of the cross-linker. This value may be negative if cross-linking produces a loss of mass.


## Output File Type #2: Percolator Input Files ##

Additionally, the PSMs are output into multiple summary files that are pre-formatted for use with [Percolator](Percolator.md). These summary files are organized by type of PSM (inter-protein cross-links, loop-links, etc.) and only contain PSMs. They are labeled as Inter, Intra, Loop, Single, and Dimer (if selected in the parameters).

  1. **SpecId** - A unique indicator label for each PSM.
  1. **Label** - Indicates target or decoy result as 1 or -1, respectively.
  1. **Scan Number** - The scan even number.
  1. **Score** - The Kojak cross correlation score. For cross-linked PSMs, this is the summed Score of each peptide.
  1. **dScore** - The delta score, or difference between the Score of this PSM and the Score of the next best PSM.
  1. **NormRank** - Inter- and Intra-protein cross-links only. The sum of the ranks of each peptide from the first pass of the analysis.
  1. **PPScoreDiff** - For Inter- and Intra-protein cross-links only. The difference between the Score and the contribution from the best scoring peptide of a cross-link, i.e. the Score of the lower scoring peptide.
  1. **Charge** - The observed precursor charge for the spectrum.
  1. **Mass** - The observed precursor mass for the spectrum.
  1. **PPM** - The mass error between Mass and PSM, in parts per million.
  1. **Len** - Single peptides and Loop-links only. The length of the peptide sequence.
  1. **LenShort** - Dimers, Inter- and Intra-protein cross-links. The length of the shorter peptide.
  1. **LenLong** - Dimers, Inter- and Intra-protein cross-links. The length of the longer peptide.
  1. **LenSum** - Dimers, Inter- and Intra-protein cross-links. The sum of the length of both peptides.
  1. **Peptide** - The peptide sequence of the PSM. When looking at Dimers and cross-links, the two peptides are represented as a single string of characters, separated by a plus (+) or two dashes (--), respectively. After each peptide sequence, the relative site of linking is provided in parenthesis. For Loop-links, both sites of linking are provided in a single set of parenthesis. Null flanking characters (-. and .-) bookend the peptide string to satisfy Percolator format rules.
  1. **Proteins** - The protein(s), tab-delimited, that contain the peptide(s) of the PSM.