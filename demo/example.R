library(QQperm)

# load igm data example
data(example.data)

#set number of permutations to generate NULL distribution of P-values
n.permutations <- 1000 #too low for real analysis, default value is 1000.

#caclualte expected and observed distributions of P-values using igm data
Ps <- igm.get.pvalues(example.data$data,example.data$is.case,n.permutations)

#write ourput to pdf file only if the R is not running in interactive mode.
if (!interactive()) {
  pdf("QQ_output.pdf")
}

#do qq plot
qqplot(Ps$perm, Ps$observed)

#estimate inflation factor lambda and plot the result
lambda <-estlambda2(Ps$observed,Ps$perm, plot = TRUE, adjust.xy = TRUE)

#write ourput to pdf file only if the R is not running in interactive mode.
if (!interactive()) {
  dev.off()
}

