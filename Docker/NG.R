#!/usr/bin/Rscript

library(ggplot2)
args<-commandArgs(TRUE)
file_input <- read.table(args[1], quote="\"", comment.char="")
NGX_test <- file_input$V1
NA_vector <- rep(NA, 100-(length(NGX_test)))
NGX_scaffoldlength <- c(NGX_test, NA_vector)
log_scaled <- log(NGX_scaffoldlength)
NGX_threshold <- c(1:100)
xdata <- NGX_threshold
NGXdata <- data.frame(log_scaled, xdata)
p = ggplot() +
geom_line(data = NGXdata, aes(x = xdata, y = log_scaled)) +
ylab("Log scaled scaffold_length") +
xlab("NG(X)") +
scale_x_continuous(breaks = seq(0, 100, by = 10)) +
geom_vline(xintercept = 50)
ggsave("NG_plot.png", plot = p)
