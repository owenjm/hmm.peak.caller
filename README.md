# hmm.peak.caller
Simple HMM-based peak caller for broad domains.

## Required packages
Requires the [RHmm](https://r-forge.r-project.org/R/?group_id=85) R package.

## Usage

```
Rscript hmm.peakcaller.R [bedgraph file]
```

Run with default settings to call three states on binding data in bedgraph format and the script will select the state with the highest median as peaks.  Works well with data in broad domains (e.g. Lamin-B, some chromatin-modifying proteins).

```
Options:                                                                        
  Number of states to fit HMM:                                                  
  --nStates=3                                                                   
                                                                                
  Number of random starts for HMM fitting:                                      
  --rhmm.iter=200                                                               
                                                                                
  Number of iterations to run per random start:                                 
  --rhmm.iter.init=10                                                           
                                                                                
  Number of cores to use for HMM fitting:                                       
  --mc.cores=8                                                                  
                                                                                
  If specified, previously saved multicore seeds will be loaded from this file: 
  --load.mcore.seeds=                                                           
                                                                                
  Chromosomes to fit HMM model to (separate by commas, no spaces):              
  --chr.model=all                                                               
```

Run with --help to see all available options.
