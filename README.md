# DIA-Umpire-Maven
DIA-Umpire, built with maven

## Refs

[authors' website](http://diaumpire.sourceforge.net/)

> DIA-Umpire is an open source Java program for computational analysis of data independent acquisition (DIA) mass spectrometry-based proteomics data.
> It enables untargeted peptide and protein identification and quantitation using DIA data, and also incorporates targeted extraction to reduce the number of cases of missing quantitation.

> As a result of recent improvements in mass spectrometry (MS),
> there is increased interest in data-independent acquisition
> (DIA) strategies in which all peptides are systematically
> fragmented using wide mass-isolation windows (‘multiplex
> fragmentation’)

- muliplexed ms/ms for improved data-independent acquisition abstract

> In mass spectrometry–based proteomics, data-independent
> acquisition (DIA) strategies can acquire a single data set
> useful for both identification and quantification of detectable
> peptides in a complex mixture. However, DIA data are noisy
> owing to a typical five- to tenfold reduction in precursor
> selectivity compared to data obtained with data-dependent
> acquisition or selected reaction monitoring. We demonstrate
> a multiplexing strategy, MSX, for DIA analysis that increases
> precursor selectivity fivefold.

- authors' 2015 paper abstract

> As a result of recent improvements in mass spectrometry (MS),
> there is increased interest in data-independent acquisition
> (DIA) strategies in which all peptides are systematically
> fragmented using wide mass-isolation windows (‘multiplex
> fragmentation’). DIA-Umpire (http://diaumpire.sourceforge.
> net/), a comprehensive computational workflow and
> open-source software for DIA data, detects precursor and
> fragment chromatographic features and assembles them into
> pseudo–tandem MS spectra. These spectra can be identified
> with conventional database-searching and protein-inference
> tools, allowing sensitive, untargeted analysis of DIA data
> without the need for a spectral library. Quantification is
> done with both precursor- and fragment-ion intensities.
> Furthermore, DIA-Umpire enables targeted extraction of
> quantitative information based on peptides initially
> identified in only a subset of the samples, resulting in
> more consistent quantification across multiple samples.
> We demonstrated the performance of the method with
> control samples of varying complexity and publicly available
> glycoproteomics and affinity purification–MS data.

- authors' 2016 paper abstract

> We describe an improved version of the data-independent acquisition (DIA) computational anal-
> ysis tool DIA-Umpire, and show that it enables highly sensitive, untargeted, and direct (spectral
> library-free) analysis of DIA data obtained using the Orbitrap family of mass spectrometers.
> DIA-Umpire v2 implements an improved feature detection algorithm with two additional filters
> based on the isotope pattern and fractional peptide mass analysis. The targeted re-extraction
> step of DIA-Umpire is updated with an improved scoring function and a more robust, semi-
> parametric mixture modeling of the resulting scores for computing posterior probabilities of
> correct peptide identification in a targeted setting. Using two publicly available Q Exactive
> DIA datasets generated using HEK-293 cells and human liver microtissues, we demonstrate
> that DIA-Umpire can identify similar number of peptide ions, but with better identification
> reproducibility between replicates and samples, as with conventional data-dependent acquisi-
> tion. We further demonstrate the utility of DIA-Umpire using a series of Orbitrap Fusion DIA
> experiments with HeLa cell lysates profiled using conventional data-dependent acquisition and
> using DIA with different isolation window widths.

## Concepts

### DIA Strategies

- multiplex fragmentation(all peptides are systematically fragmented using wide mass-isolation windows)

### DIA-Umpire Computations

- DDA(data-dependent acquisition), DIA(data-independent acquisition)

- MS1, MS2 data detects monoisotopic masses and elution profile shapes.

- here featuire is about the chromatographic feature

- dia-spectra(MS1, MS2 data) after feature detection -> possible precursor and fragment ion features

- assembles them into pseudo–tandem MS spectra

- pseudo MS/MS spectra (from MS1 features grouped with fragments) for untargeted MS/MS database search to identify peptides and proteins

- protein and peptide identifications, confidently identified peptides

- These spectra can be identified with conventional database-searching and protein-inference tools,
allowing sensitive, untargeted analysis of DIA data without the need for a spectral library

- the analysis of MS/MS database search is not provided by DIA-Umpire, Trans Proteomics Pipeline (TPP) for this purpose.

- untargeted MS/MS search

- dia-umpire targeted re-extraction

- building an internal spectral library.

- All IDs from either untargeted MS/MS database search or targeted re-extraction are linked to the corresponding precursor-fragment groups which carry quantitative information in the form of precursor and fragment ion intensities. This quantitative information is stored at different levels (fragment → peptide → protein) and is reported in Step D.
