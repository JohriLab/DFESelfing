# a script to calculate pi using m1 mutations
# in slim output files

#m1_files come from grepping all mutations with "m1" from files. Command was:
#for x in *DFE[[:digit:]]; do cd $x; for y in *.txt; do grep "m1" $y > m1_$y; done; cd ..; done

#calculation will be p = last_column/100
#total_neu_mutations = 187500
#sum(2*p*(1-p))/total_neu_mutations

rm(list=ls())

total_neu_mutations <- 187500

#need to read in all m1_output.txt files in a way that preserves selfing%_DFE#_replicate structure
#get pi for each output
#average and std dev for each DFE and make table like I did for gamma and beta

dirs <- dir(path = "/nas/longleaf/home/adaigle/adaigle/johri_elegans/gammaDFE", pattern = "DFE")
dirs <- list.dirs(path = "/nas/longleaf/home/adaigle/adaigle/johri_elegans/gammaDFE")
dirs <- dirs[ grepl("eqm", dirs) ]

# has code to load in files with nice names, but i realized I didn't need that 
#makes a vector of avg pi for each experiment 
#outputs are averaged and std deviation put into table
summary_table <- data.frame()
for(x in dirs) {
    replicate_pis <- c()
    for(y in list.files(path = x, pattern = "m1")) {
        #assign( 
        #paste(sub(
        #    sub("[^/]+$", "", x),
        #   "", x), y, sep=""), read.table(paste(x, "/", y, sep = "")))
        z <- read.table(paste(x, "/", y, sep = ""))
        replicate_pis <- c(replicate_pis, sum(2*(z$V9/100) *(1-(z$V9/100)))/total_neu_mutations)
        #avg <- mean z$V10
    }
    summary_table <- rbind(summary_table, c(sub(sub("[^/]+$", "", x),"", x),
        mean(replicate_pis), sd(replicate_pis)))
}

row.names(summary_table) <- summary_table[[1]]
summary_table <- summary_table[-1]
colnames(summary_table) <- c("pi", "pi_sd")
#relative biases due to selfing/dfe
#E(pi) is theta in a neutral model
theta <- 4*5000*3.3e-9*100
#aka B, background selection
summary_table$B <- as.numeric(summary_table$pi)/(theta)
#empircal Ne = B*Ne
summary_table$empirical_Ne <- summary_table$B*5000
write.csv(summary_table, file="/nas/longleaf/home/adaigle/DFESelfing/pi_summary.csv")



