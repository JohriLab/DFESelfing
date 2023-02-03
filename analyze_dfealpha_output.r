path_to_outputs <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output/")

stats <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
colnames(stats) <- "stats"
DFE_list <- c("DFE1", "DFE2", "DFE3")
replicates <- c("output1", "output2","output3","output4","output5")


for(x in DFE_list) {
    for(y in replicates) {
        data <- readLines(
            paste("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output/", 
            x, y, ".txt_selected/est_dfe.out", sep = ""))
        vector_data <- unlist(strsplit(data, " "))
        df <- data.frame(matrix(vector_data, ncol = 2, byrow = TRUE))
        colnames(df) <- c("stats", paste(x,y,sep=""))
        stats <- merge(stats,df)
    }
}

stats <- as.data.frame(t(stats))
names(stats) <- stats[1,]
stats <- stats[-1,]
stats$NwEs <- as.numeric(stats$Es) * as.numeric(stats$Nw)
stats$gamma <- 2*stats$NwEs

means <- data.frame()
for(x in DFE_list) {
    tmpdf <- stats[grep(x, rownames(stats)),]
    tmpdf <- as.data.frame(sapply(tmpdf, FUN=as.numeric))
    means <- rbind(means, c(x, sapply(tmpdf, FUN=mean)))
}
colnames(means) <- c("Experiment", colnames(stats))
