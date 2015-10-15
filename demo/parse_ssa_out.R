source('read_csv_list.R')

out_summary <- function(file) {
    tbls <- read.csv.list(file)
    emp <- data.frame(table(factor(tbls$samples$k)))
    names(emp) <- c('k','emp')
    emp$emp <- emp$emp/sum(emp$emp)
    ps <- tbls$propensities
    ps$k <- factor(ps$k)
    ps$p <- ps$p/sum(ps$p)
    merge(ps,emp,by='k')
}

plot_out_summary <- function(os,...) {
    miny <- min(0,os$p, os$emp)
    maxy <- max(os$p, os$emp)
    ks <- as.vector(os$k)
    plot(ks, os$p, ylim=c(miny,maxy),...)
    points(ks, os$emp, pch='x',...)
}

