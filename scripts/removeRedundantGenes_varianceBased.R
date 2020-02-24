args <- commandArgs(trailingOnly= T)
#file.format <- args[1] ## Either: 1) _Roadmap or 2) and empty string
#data.name <- args[2]
#file.path <- args[3]

feat.input <- args[1]
resp.input <- args[2]
feat.output <- args[3]
resp.output <- args[4]

#file.format <- "_Roadmap" ## Alternative" empty string to refer to the TF affinities, I originally started the project based on
#file.format <- ""


#feat <- read.table(paste(file.path, "/scMTL_", data.name, "_feature", file.format, ".txt", sep =""), header=T, row.names= 1)
#response <- read.table(paste(file.path, "/scMTL_", data.name, "_response", file.format, ".txt", sep= ""), header=T, row.names= 1)

feat <- read.table(feat.input, header=T, row.names= 1)
response <- read.table(resp.input, header=T, row.names= 1)

#feat.vars <- apply(as.matrix(feat), 1, FUN= var)
exp.vars <- apply(as.matrix(response), 1, FUN= var)

#cutoff <- summary(feat.vars)[5]
## With some datasets the 5th quantile was rsulting in very few genes to survive the filter, therefore I added this min operator that controls what value to be chosen as the cutoff
cutoff <- min(1, summary(exp.vars)[5])

#low.var.idx <- which(feat.vars < cutoff)
low.var.idx <- which(exp.vars < cutoff)
#low.var.idx <- which(exp.vars < .05) # We used this criteria for the monocle data that had 306 cells. As the number of cells increase, statistically speaking the variance should decrease and this fixed cutoff does not make sense anymore. Therefore, I'm gonna enable the data-driven cutoff by uncommenting the line above, especially for bigger datasets like StemNet from Kathrin


## Commenting the discarded parts out for the snakemake purposes!
#write.csv(feat[low.var.idx, ], file= paste(file.path, "/scMTL_", data.name, "_feature_discarded_var", file.format, ".txt", sep= ""))
#write.csv(response[low.var.idx, ], file= paste(file.path, "/scMTL_", data.name, "_response_discarded_var", file.format, ".txt", sep= ""))

write.csv(feat[-low.var.idx, ], file= feat.output)
write.csv(response[-low.var.idx, ], file= resp.output)

#write.csv(feat[-low.var.idx, ], file= paste(file.path, "/scMTL_", data.name, "_feature_reduced_var", file.format, ".txt", sep= ""))
#write.csv(response[-low.var.idx, ], file= paste(file.path, "/scMTL_", data.name, "_response_reduced_var", file.format, ".txt", sep= ""))

