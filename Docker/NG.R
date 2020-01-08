#!/usr/bin/Rscript

library(ggplot2)
library (optparse)
library(reshape)
option_list <- list ( make_option (c("-f","--filelist"),
                                   help="comma separated list of files (default %default)")
                     )

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

myfilelist <- strsplit(opt$filelist, ",")
M <- unlist(myfilelist)


list1 <- list()
for(i in 1:length(M)){
    file <- read.table(M[i], quote="\"", comment.char="")
    p <- length(file$V1)
    NA_vector <- rep(NA, 100-p)
    #print (NA_vector)
    #NGX_test[i] <- list()
    h <- as.vector(file$V1)
    NGX_scaffoldlength <- c(h, NA_vector)
    #k <- as.vector(NGX_scaffoldlength)
    list1[[i]] <- NGX_scaffoldlength
}

l <- as.vector(list1)


X1 <- as.data.frame(l)
NGX_threshold <- c(1:100)
colnames(X1) <- M
xdata <- NGX_threshold
NGXdata <- data.frame(xdata, X1)
xymelt <- melt(NGXdata, id.vars = "xdata")


p <- ggplot(xymelt, aes(x = xdata, y = value, color = variable)) +
     theme_bw() +
     geom_line() +
     ylab("scaffold_length (bp)")+
     xlab("NG(X)") +
     scale_x_continuous(breaks = seq(0, 100, by = 10)) +
     geom_vline(xintercept = 50)
ggsave("NG_plot.png", plot = p)

