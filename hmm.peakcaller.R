# Quick and dirty HMM fitting to bedgraph / GFF bound genomic datasets
# Copyright © 2020, Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA

library(RHmm,quietly=T)
library(parallel)
library(tools)

in.files = vector()

usage.message = "\ndamid.hmm.peakcaller -- a quick and dirty HMM-based peak caller for genomic binding data.\nWorks best to call peaks on proteins with broad chromatin-binding domains.\n\n"

op.args = list(
  "nStates" = 3,
  "rhmm.iter" = 200,
  "rhmm.iter.init" = 10,
  "mc.cores" = 8,
  "load.mcore.seeds" = "",
	"chr.model" = "all"
)

op.notes = list(
  "nStates" = "Number of states to fit HMM",
  "rhmm.iter" = "Number of random starts for HMM fitting",
  "rhmm.iter.init" = "Number of iterations to run per random start",
  "mc.cores" = "Number of cores to use for HMM fitting",
  "load.mcore.seeds" = "If specified, previously saved multicore seeds will be loaded from this file",
	"chr.model" = "Chromosomes to fit HMM model to (separate by commas, no spaces)"
)

### save random seed for future reproducibility
test = runif(1)
my.seed  = .Random.seed
write.table(my.seed,".randomseed")

### Read CLI options
input.args = commandArgs(trailingOnly = TRUE)

help = function () {
  cat(usage.message)
  cat("Options:\n")
  for (n in names(op.args)) {
    cat(paste("  ",op.notes[[n]],":\n",sep=""))
    cat(paste("  --",n,"=",op.args[[n]],"\n\n",sep=""))
  }
  cat("\n")
  quit("no",1)
}

read.ops = function (x) {
  for (op in x) {
    if (any(grepl("^--",op))) {
      op = gsub("^--","",op)
      y = unlist(strsplit(op,"="))
  
      if (y[1] == "help") {
        help()
      }
  
      if (!is.null(op.args[[ y[1] ]])) {
        op.args[[ y[1] ]] <<- y[2]
      } else {
        cat("Error: Option",y[1],"not recognised ...\n")
        quit("no",1)
      }
    } else {
      in.files <<- c(in.files,op)
    }
  }
}

write.ops = function () {
  oldw = getOption("warn")
  options(warn = -1)
  out.df = data.frame(option="version",value=my.version)
  for (n in names(op.args)) {
    v <<- op.args[[n]]
    df.line = data.frame(
      option=n,
      value=v
    )
    out.df = rbind(out.df, df.line)
  }
  write.table(out.df,"input.args.txt",row.names=F)
  options(warn = oldw)
}

read.ops(input.args)

condense.gff = function (input.df) {
	total = nrow(input.df)
	merged = list()
	
	chr.last = ""
	state.start = NA
	state.end = NA
	state.current = NA
	rcount = 1
		
	for (i in c(1:total)) {
	  if (i%%1000 == 0) {cat(paste("Processing row",i,"of",total,"             \r"))}
	  
	  chr = .subset2(input.df,1)[i]
	  start = .subset2(input.df,2)[i]
	  end = .subset2(input.df,3)[i]
	  state = .subset2(input.df,4)[i]
	  
	  if (chr == chr.last) {
		if (state == state.current) {
		  state.end = end
		} else {
		  # save state block
		  merged[[rcount]] = data.frame(chr=chr, start=state.start, end=state.end, state=state.current)
		  rcount = rcount+1
		  
		  state.start = start
		  state.end = end
		  state.current = state
		}
	  } else {
		# new chromosome
		state.start = start
		state.end = end
		state.current = state
	  }
	  
	  chr.last = chr
	}
	merged[[rcount]] = data.frame(chr=chr, start=state.start, end=state.end, state=state.current)
	
	merged.df = do.call('rbind',merged)
	return(merged.df)
}

convert.hex.rgb = function (x) {
  hex.r = substr(x,2,3)
  hex.g = substr(x,4,5)
  hex.b = substr(x,6,7)
  dec.r = strtoi(hex.r,16L)
  dec.g = strtoi(hex.g,16L)
  dec.b = strtoi(hex.b,16L)
  out = paste(dec.r,",",dec.g,",",dec.b,sep="")
  return(out)
}

make.states.bed = function (datf) {
  # expect "states gff" input: chr start end state
  datf$zero = 0
  datf$strand = "."
  datf$cols = sapply( datf$state, function(x) convert.hex.rgb(gg_colour_hue(nStates)[x]) )
  
  states.bed = data.frame(datf$chr,datf$start,datf$end,datf$state,datf$zero,datf$strand,datf$start,datf$end,datf$cols)
  return(states.bed)
}

fit.hmm.mc = function (data = data.na, cores = mc.cores, chrs = NULL) {
  # multi.core HMM to speed up model fitting
  
  # Allow restriction of training dataset to individual chromosome or chromosomes
  # Reducing the training set will decrease time to fit model (but use with some caution)
  
  model.input = data.frame()
  if (chrs[1] == "all") {
    model.input = data
  } else {
    for (c in chrs) {
      model.input = rbind(model.input,data[data$chr == c,])
    }
  }
  
  cat(paste("  Fitting",nStates,"states ...\n"))
  
  # We need different random seeds for each thread, but want this operation to be reproducible
  seeds = sample(0:2147483647,mc.cores,replace=F)
  if (op.args[["load.mcore.seeds"]] != "") {
    if (file.exists(op.args[["load.mcore.seeds"]])) {
      # load previously saved random seeds
      cat("  Loading seeds file ...\n")
      seeds = read.table(paste("../",op.args[["load.mcore.seeds"]],sep=""))[[1]]
    } else {
      cat("  Error, cannot read seeds file.\n")
    }
  }
  write.table(seeds,"multicore.random.seeds",sep="\n",quote=F,row.names=F,col.names=F)
  
  # Fit multi-threaded
  mc.model = mclapply(seeds, function (x) {
      set.seed(x)
      HMMFit(
      data.matrix(model.input[,4:ncol(model.input)]),
      nStates=nStates,
      control=list(
                   verbose=2,
                   nInit=as.integer(rhmm.n.iter/length(seeds)),
                   nIterInit=rhmm.n.iter.init
                   )
      )},
      mc.cores=cores
  )
  
  hmm.best = vector()
  for (i in 1:length(mc.model)) {
    if (is.finite(mc.model[[i]]$LLH)) {
      # We should probably throw a warning if NAs are introduced into LLH values, but for now we're just ignoring them ...
      hmm.best=c(hmm.best,mc.model[[i]]$LLH)
    }
  }
  
  cat("Model fitted.\n")
  print(hmm.best)
  print( grep(max(hmm.best),hmm.best) )
  
  return( mc.model[[ grep(max(hmm.best),hmm.best) ]] )
}

read.gff = function (x,name="score") {
  fn.ext = file_ext(x)
  
  if (grepl("gff",ignore.case=T,fn.ext)) {
	temp.data <- read.table(x,row.names=NULL)
  	if (ncol(temp.data) > 5) {
  	  # GFF
  	  trim.data = temp.data[,c(1,4,5,6)]
  	} else {
  		cat("Error: file does not appear to be in GFF format\n\n")
  		quit("no",1)
  	}
  } else if (grepl("bed",ignore.case=T,fn.ext)) {
  	temp.data = read.table(x,row.names=NULL,skip=1)
  	if (ncol(temp.data) == 4) {
  		# bedgraph
  		trim.data = temp.data
  	} else {
  		cat("Error: file does not appear to be in bedGraph format\n\n")
  		quit("no",1)
  	}
  } else {
  	cat("Error: input file does not appear to be in bedGraph or GFF format ...\n\n")
  	quit("no",1)
  }
  
  names(trim.data) = c("chr","start","end",name)
  trim.data$chr = gsub("^chr","",trim.data$chr,perl=T)
  
  return(trim.data)
}

gg_colour_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

if (length(in.files) == 0) {
  help()
}

# set globals
filename = in.files[1]

nStates = as.integer(op.args[["nStates"]]) # number of HMM states to fit
rhmm.n.iter = as.integer(op.args[["rhmm.iter"]])
rhmm.n.iter.init = as.integer(op.args[["rhmm.iter.init"]])
mc.cores = as.integer(op.args[["mc.cores"]])
chr.model = strsplit(op.args[["chr.model"]],",")[[1]]

cat(paste("Reading",filename,"to fit",nStates,"states ...\n"))

# read data
data = read.gff(filename)
data[data == 0] = NA
data.bedgraph = na.exclude(data)
scores = data.bedgraph[,4]

# Fit HMM
cat("Fitting HMM ...\n")
hmm = fit.hmm.mc(data = data.bedgraph, cores = mc.cores, chrs = chr.model)

cat("Generating viterbi path ...\n")
vit <- viterbi(hmm,scores)

# Find peaks
out.bedgraph.all <- data.frame(data.bedgraph,vit$states)

# peak.index is defined by the state with the highest mean.  This will cause issues if large numbers of states are fitted.
state.means = vector();
for (i in c(1:nStates)) {
  state.means[i] = mean(out.bedgraph.all$score[out.bedgraph.all$vit.states == i])
}
peak.index = grep(max(state.means),state.means)

# generate files at both GATC resolution and condensed to homogenous regions
out.bedgraph = data.frame(data.bedgraph[,1:3],vit$states)
out.cond.bedgraph = condense.gff(out.bedgraph)

out.peaks = out.cond.bedgraph[out.cond.bedgraph$state==peak.index,]
out.peaks.decond = out.bedgraph[out.bedgraph$vit.states==peak.index,]

cat("Writing data ...\n")
write.table(file=paste(filename,"peaks","gff",sep="."),data.frame(out.peaks$chr,".",".",out.peaks$start,out.peaks$end,1,".",".","."),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file=paste(filename,"peaks.decond","gff",sep="."),data.frame(out.peaks.decond$chr,".",".",out.peaks.decond$start,out.peaks.decond$end,1,".",".","."),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

out.cols.bed = make.states.bed(out.cond.bedgraph)
write('track name="" description="" visibility=2 itemRgb="On"', file=paste(filename,".rhmm_",nStates,"_states.cols.bed",sep=""))
write.table(out.cols.bed,file=paste(filename,".rhmm_",nStates,"_states.cols.bed",sep=""),sep="\t",append=T,row.names=FALSE,col.names=FALSE,quote=FALSE)
  
