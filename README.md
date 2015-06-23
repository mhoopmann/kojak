# Kojak

A database search algorithm for solving cross-linked peptide mass spectra.

## License

Copyright 2014-2015, Michael R. Hoopmann, Institute for Systems Biology

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
  
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## Kojak Release Notes

#### Version 1.3.6 - June 17 2015
* Fixed bug where protein site of linkage was offset by 1.
* Fixed bug where heterobifunctional cross-links were duplicated in the results.
* Fixed bug when reporting negative cross-link masses.
* Fixed buffer overrun when pairing peptide mass to spectrum mass range.
* Added toggle for searching differential modifications on cross-linked peptides.
* Added toggle for searching mono-linked sites on cross-linked peptides.
* Improved precursor refinement algorithm.
* Extended prefer_precursor_pred to supplement precursor info with Kojak determined precursors.
* Extended diagnostics to support multiple spectra.

#### Version 1.3.5 - April 15 2015
* Fixed bug in setting modifications to protein C-terminus.
* Fixed Inter/Intra Percolator-formatted results for real this time.
* Added amine+ reactive group for cross-linking (lysines and n-terminus, plus alternate sites of serine, threonine, and tyrosine).

#### Version 1.3.4
* Added cross-linking of acyl groups.
* Fixed bug in labeling of Inter and Intra cross-linked Percolator files.

#### Version 1.3.3
* Percolator results now divided into multiple files based on PSM type.

#### Version 1.3.2
* Improved Kojak Xcorr algorithm to use less memory
* Faster Xcorr implementation over Version 1.0
* Support for multi-threaded systems

#### Version 1.0
* Initial Release
