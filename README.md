# hmm.peak.caller
Simple HMM-based peak caller for broad domains.

## Required packages
Requires the [RHmm](https://r-forge.r-project.org/R/?group_id=85) R package.

## Usage
Run with default settings to call three states on binding data in bedgraph format and the script will select the state with the highest median as peaks.  Works well with data in broad domains (e.g. Lamin-B, some chromatin-modifying proteins).

Run with --help to see all available options.
