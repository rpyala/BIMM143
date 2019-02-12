#' ---
#' title: "Lecture 5 Lab Portion"
#' author: "*Reshma Pyala*"
#' date: "*January 24th, 2019*"
#' output: pdf_document
#' ---

#' This is some text and we can **bold**, *italic*, and `code` words so that it may appear like so in our document. 


# Class 05 R graphics intro

#My first boxplot
x <- rnorm(1000,0)
boxplot(x)

summary(x)
hist(x)

boxplot(x, horizontal = TRUE)


#Lab portion - Line Plot (2A)
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)

plot(weight[,1], weight[,2], pch=2, cex=1, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", typ="o", main="Baby weight with age", lty=1)

#Lab Portion - Bar Plot (2B) - COME BACK TO THIS PORTION
#use the "\t" to get rid of the error and changing the header will likely get rid of the error message
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", sep="\t", header=TRUE)
par()$mar
par(mar=c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, names.arg = mouse$Feature, horiz=TRUE, ylab = "", main = "Number of features in the mouse GRCm38 genome", las=1, xlim = c(0, 80000), col=rainbow(11))
#las give helps with labeling, las=1 means always horizontal labels, par()$mar helps find the margins and dimensions of the graph

#Lab portion - Histograms
c(rnorm(10000), rnorm(10000)+4, breaks=10)
x <- c(rnorm(10000), rnorm(10000)+4)
hist(x, breaks=80, xlim = c(-4, 8), ylim= c(0, 900))

#Lab portion - Using color in plots
read.table("bimm143_05_rstats/male_female_counts.txt", sep="\t", header=TRUE)
mf <- read.table("bimm143_05_rstats/male_female_counts.txt", sep="\t", header=TRUE)
barplot(mf$Count, names.arg = mf$Sample, col=rainbow(nrow(mf)), las=2, ylab="Counts")

barplot(mf$Count, names.arg = mf$Sample, col=c("blue3", "red3"), las=2, ylab="Counts")

#next portion of color - colouring by value
read.table("bimm143_05_rstats/up_down_expression.txt", "\t", header=TRUE)
genes <- read.table("bimm143_05_rstats/up_down_expression.txt", "\t", header=TRUE)
nrow(genes)
#Has 5196 genes

#How many up, down, and all around?
table(genes$State)
#down = 72, unchanging = 4997, up = 127
plot(genes$Condition1, genes$Condition2, col=genes$State, xlab="Expression Condition1", ylab="Expression Condition2")
palette()
levels(genes$State)
palette(c("blue", "grey", "red"))
plot(genes$Condition1, genes$Condition2, col=genes$State, xlab="Expression Condition 1", ylab="Expression Condition 2", pch=16)
