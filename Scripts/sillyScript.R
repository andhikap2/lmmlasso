#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

## program...
df = read.table(args[1], header=TRUE)
num_vars = which(sapply(df, class)=="numeric")
df_out = df[ ,num_vars]
write.table(df_out, file=args[2], row.names=FALSE)