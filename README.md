# Automated-Single-Molecule-Colocalization-Analysis
Analysis of Single molecules binding to DNA fragments using data from TIRF microscopy- Biggins lab

Automated-Single-Molecule-Colocalization-Analysis (ASMCA) repository contains the MATLAB (R2019b) codes used to generate the data related to TIRF microscopy experiments presented in Popchock *et al.* from the Biggins Lab. 

The aim of this pipeline is to identify binding events between DNA fragments immobilized on coverslips incubated with a cellular extract containing fluorescently tagged proteins and quantify their duration using time-lapse TIRF microscopy.
In short,  the image dataset was drift-corrected using either fast Fourier Transform cross-correlations or translation affine transformation depending on the severity of the drift. DNA spots were identified after binarizing the DNA signal using global background value as threshold, followed by size and time-persistency filtering. Mean values of z-normalized fluorescent markers intensities were measured at each DNA spot at each time frame, and local background was subtracted. Z-normalized traces were then binarized to ON/OFF pulses by applying a channel-specific, manually adjusted threshold value unique to all traces in a given image set. Pulses onsets, durations and overlaps between channels were then derived. Pulses in ON state at the beginning or the end of the recording were censored. See codes and manuscipt for detail. Input data is a three-channel time-series generated from a Nikon TIRF microscope with settings and specifications described in the manuscript. This pipeline has not been tested on other datasets coming from different platforms.

*tirfregister.m* corrects the image dataset for the image drifting during timelapse acquisition.

*ascma.m* analyzes the colocalization events occuring at the DNA fragments. For image dataset import, it relies on the bio-format package for MATLAB *bfmatlab* that can be freely downloaded at:

https://www.openmicroscopy.org/bio-formats/downloads/

*FastPeakFind.zip* contains *FastPeakFind.m*, a function created by Adi Natan that is used to find local maxima in noisy images. It is called by the main *ascma* function to find centroids of DNA spots.

*progressText.m* is a helper function developped by the Danuser lab to show progress of a loop as text on the screen. It is employed by both *tirfregister()* and *ascma()*.

*ExampleRaw.tif* is a cropped, shortened dataset derived from an original .nd2 image dataset that can be used for demo purposes.

```
ISOR = tirfregister('ExampleRaw.tif',0);
```

*ExampleReg.ome.tiff* is the drift-corrected version of *ExampleRaw.tif* that can be used for demo:

```
[data, SK, ISOR, StatsPix, survData] = asmca('ExampleReg.ome.tiff', 5, 0.4, 0.4);
```

