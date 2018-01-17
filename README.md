# Pipeline Overview

## 1. Setup and configuration
A very straightforward file called "config" (suffix .m, .py, or .R as appropriate) should be placed in the code/matlab, code/python, and code/R. The config file needs to specify three variables: data_dir, cache_dir, and plot_save_dir. Each should be a path pointing to the data directory, a location avaialable for caching, and the directory to save plots.
## 2. Data cleaning
The fundamental data used in the MATLAB portions are structs with a structure specific to the data collection. Not all of it is needed, so a good deal of the metainformation is discarded in creating generalized holders for interacting with the underlying data (see DataHandler).
### 2.1 Eliminating electrode outliers
In the initial creation of a DataHandler, individual electrodes are (optionally) eliminated from the data set due to their extreme dissimilarity with the other electrodes. Mean cosine simmilarity with each other electrode is used as the simmilarity measure. Usually only a couple additional electrodes are excluded in this process. Note cosine simmilarity is calculated only with LFPs, but spike data associated with the "bad" electrodes are also eliminated.
### 2.2 Sorting spikes
The initial structuring of the data is trial-based, so there are separate fields for a trial and intertrial period. Even though spikes have global timestamps (i.e. not relative to a given trial or intertrial period), spikes tend not to be sorted across the trial sections. To facilitate analysis and alleviate headaches, spikes are sorted when they are combined in a single cell array in the DataHandler.
### 2.3 Subsetting data from existing DataHandlers
Lastly, there are many readily apparent anoalies in the LFP data (very large and highly rhythmic spikes). The data was therefore carefully selected to avoid the aberrant behavior. The exact times to split were determined via visual inspection.
## 3. Analysis
The DataHandler class provides a couple methods of subsetting data based on trial conditions. Using these, it is possible to perform analysis on a variety of trial conditions and lengths of data. For example, LFP data can easily be selected from between noise onset and shape stimulus events and then passed to various spectral analysis methods (provided by Chronux). Using the electrode mapping information kept in each DataHandler, it is relatively easy to calculate LFP coherences and bin them based on the distance between the two electrodes. 
## 4. Plotting

# Code Summary

## 1. DataHandler
### 1.1 Fields
#### Data
trials
spikes
lfps
#### Trial
fixate
noise
shape
saccade

### 1.2 Methods
Constructor
select
getDataSlices

## 2. Analysis
## 3. Plotting
## 4. Utility methods
### 4.1 General data manipulation
### 4.2 Conversion between indices
