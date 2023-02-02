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
data <- readLines("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output/DFE1output1.txt_selected/est_dfe.out")
vector_data <- unlist(strsplit(data, " "))


# create a vector of the data

# convert the vector to a data frame
df <- data.frame(matrix(vector_data, ncol = 2, byrow = TRUE))

merge(stats,df)
