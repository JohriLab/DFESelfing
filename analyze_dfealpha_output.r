rm(list=ls())
library(RColorBrewer)

#create a dataframe with means for replicates of each experiment
#also finds sd for gamma and beta from the replicates
find_means <- function(path_to_outputs, DFE_list, replicates) {
    stats <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
    colnames(stats) <- "stats"

    for(x in DFE_list) {
        for(y in replicates) {
            data <- readLines(
                paste(path_to_outputs, 
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
    sdbeta <- c()
    sdgamma <- c()
    for(x in DFE_list) {
        tmpdf <- stats[grep(x, rownames(stats)),]
        tmpdf <- as.data.frame(sapply(tmpdf, FUN=as.numeric))
        means <- rbind(means, c(x, sapply(tmpdf, FUN=mean)))
        sdbeta <- c(sdbeta, sd(as.numeric(tmpdf$b)))
        sdgamma <- c(sdgamma, sd(as.numeric(tmpdf$gamma)))

    }
    colnames(means) <- c("Experiment", colnames(stats))
    means$betasd <- sdbeta
    means$gammasd <- sdgamma
    return(means)
}

DFE_list <- c("DFE1", "DFE2", "DFE3")
replicates <- c("output1", "output2","output3","output4","output5")
selfing0 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output/")
selfing50 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output_50/")
selfing99 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output_100/")

#run function to create dataframes for each selfing percentage
means0 <- find_means(selfing0, DFE_list, replicates)
means50 <- find_means(selfing50, DFE_list, replicates)
means99 <- find_means(selfing99, DFE_list, replicates)

#betas and gammas used to generate DFE1,2,3
truebeta <- c(0.9,0.5,0.3)
truegamma <- c(5,50,1000)

#makes a summary table showing accuracy of beta and gamma predictions
gammabetatable <- data.frame(truebeta, means0$b, means0$betasd, means50$b, means50$betasd, 
means99$b, means99$betasd, truegamma, -as.numeric(means0$gamma), means0$gammasd, 
-as.numeric(means50$gamma), means50$gammasd, -as.numeric(means99$gamma), means99$gammasd)

#makes figures showing beta and gamma predictions vs truth
#will not be used because of the large difference, table is better
test <- t(matrix(c(truebeta,as.numeric(means0$b),as.numeric(means50$b), as.numeric(means99$b)), nr=3))
colnames(test) <- c("DFE1", "DFE2", "DFE3")
barplot(test, beside=T, ylab = "beta parameter",
    col=c("black", "red","blue","green"))
legend("topleft", c("truth","selfing0","selfing50", "selfing99"), pch=15, 
       col=c("black", "red","blue","green"), 
       bty="n")

threebeta <- t(matrix(c(truebeta,as.numeric(means0$b),as.numeric(means50$b)), nr=3))
colnames(threebeta) <- c("DFE1", "DFE2", "DFE3")

barplot(threebeta, beside=T, ylab = "beta parameter",
    col=c("black", "red","blue"))
legend("topleft", c("truth","selfing0","selfing50"), pch=15, 
       col=c("black", "red","blue"), 
       bty="n")

barplot(t(matrix(c(truegamma,-as.numeric(means0$gamma),-as.numeric(means50$gamma)), nr=3)), beside=T,
    col=c("black",col=c("black", "red","blue")))
legend("topleft", c("truth","selfing0","selfing50"), pch=15, 
       col=c("black", "red","blue"), 
       bty="n")

barplot(t(matrix(c(truegamma,-as.numeric(means0$gamma),-as.numeric(means50$gamma), -as.numeric(means99$gamma)), nr=3)), beside=T,
    col=c("black", "red","blue","green"))
legend("topleft", c("truth","selfing0","selfing50","selfing99"), pch=15, 
       col=c("black", "red","blue","green"), 
       bty="n")

#------------------------------------
#This is to get the proportion of f0, f1, f2 and f3 in case of diploid and haploid:

DFE_proportions <- function(Nw,meanGamma,meanS,beta) {
# this code adapted from get_dfe_class_proportions
s_shape <- beta
s_rate <- s_shape/abs(meanS)

print("gamma distribution with shape parameter, alpha:")
print(s_shape)
print("gamma distribution with rate parameter, beta:")
print(s_rate)
#these statements refer to the wikipedia way of describing gamma
#https://en.wikipedia.org/wiki/Gamma_distribution
#not how r categorizes them!

print ("DFE classes in terms of Nws: ")
x1 <- 1/(Nw)
x10 <- 10/(Nw)
x100 <- 100/(Nw)
f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
f3 <- 1.0 - f2 - f1 - f0
print(paste(f0, f1, f2, f3, sep=", "))

print ("DFE classes in terms of 2Nws: ")
x1 <- 1/(2*Nw)
x10 <- 10/(2*Nw)
x100 <- 100/(2*Nw)
f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
f3 <- 1.0 - f2 - f1 - f0
print(paste(f0, f1, f2, f3, sep=", "))
return(c(list(f0), list(f1), list(f2), list(f3)))
}

#not quite sure what the control Nw should be. Putting 100 for now
Nw <- c(as.numeric(means0$Nw), as.numeric(means50$Nw), as.numeric(means99$Nw))
#meanS <- as.numeric(args[2]) # shape = mean * rate
meanGamma <- c(-as.numeric(means0$gamma), -as.numeric(means50$gamma), -as.numeric(means99$gamma))
meanS <- meanGamma/(2.0*Nw)
beta <- c(as.numeric(means0$b), as.numeric(means50$b), as.numeric(means99$b)) #shape parameter

#using 2NwS for gamma for now, wasn't quite sure what to do. 
#vector of lists, f0, f1,f2, f4 (discrete mutation classes)
#within each class is selfing levels(DFE1-3)
#so 0selfing DFE1,2,3 ; 50 selfing DFE1-3 etc.
discrete_classes<-DFE_proportions(Nw,meanGamma,meanS,beta)

B <- (c(
0.94919771, 0.93074411, 0.89399308, 0.70231493, 0.67842912, 0.61559224, 0.16699000, 0.04243927, 0.04826596
))
#also get vectors for three true DFEs
## DFE1 true discrete
#followed by three adjustments based on pi calculations for DFE
DFE1Nw <- 100
DFE1Gamma <- 5
DFE1meanS <- DFE1Gamma/(2*DFE1Nw)
DFE1beta <- 0.9
DFE1_truth <- DFE_proportions(DFE1Nw,DFE1Gamma,DFE1meanS,DFE1beta)

DFE1meanS_S0 <- (5*0.94919771)/(2*DFE1Nw)
DFE1_truth_S0 <- DFE_proportions(DFE1Nw,DFE1Gamma,DFE1meanS_S0,DFE1beta)

DFE1meanS_S50 <- (5*0.70231493)/(2*DFE1Nw)
DFE1_truth_S50 <- DFE_proportions(DFE1Nw,DFE1Gamma,DFE1meanS_S50,DFE1beta)

DFE1meanS_S99 <- (5*0.16699000)/(2*DFE1Nw)
DFE1_truth_S99 <- DFE_proportions(DFE1Nw,DFE1Gamma,DFE1meanS_S99,DFE1beta)


## DFE2 true discrete
DFE2Nw <- 100
DFE2Gamma <- 50
DFE2meanS <- DFE2Gamma/(2*DFE2Nw)
DFE2beta <- 0.5
DFE2_truth <- DFE_proportions(DFE2Nw,DFE2Gamma,DFE2meanS,DFE2beta)

DFE2meanS_S0 <- (50*0.93074411)/(2*DFE2Nw)
DFE2_truth_S0 <- DFE_proportions(DFE2Nw,DFE2Gamma,DFE2meanS_S0,DFE2beta)
DFE2meanS_S50 <- (50*0.67842912)/(2*DFE2Nw)
DFE2_truth_S50 <- DFE_proportions(DFE2Nw,DFE2Gamma,DFE2meanS_S50,DFE2beta)
DFE2meanS_S99 <- (50*0.04243927)/(2*DFE2Nw)
DFE2_truth_S99 <- DFE_proportions(DFE2Nw,DFE2Gamma,DFE2meanS_S99,DFE2beta)

## DFE3 true discrete
DFE3Nw <- 100
DFE3Gamma <- 1000
DFE3meanS <- DFE3Gamma/(2*DFE3Nw)
DFE3beta <- 0.3
DFE3_truth <- DFE_proportions(DFE3Nw,DFE3Gamma,DFE3meanS,DFE3beta)

DFE3meanS_S0 <- (1000*0.93074411)/(2*DFE3Nw)
DFE3_truth_S0 <- DFE_proportions(DFE3Nw,DFE3Gamma,DFE3meanS_S0,DFE3beta)
DFE3meanS_S50 <- (1000*0.67842912)/(2*DFE3Nw)
DFE3_truth_S50 <- DFE_proportions(DFE3Nw,DFE3Gamma,DFE3meanS_S50,DFE3beta)
DFE3meanS_S99 <- (1000*0.04243927)/(2*DFE3Nw)
DFE3_truth_S99 <- DFE_proportions(DFE3Nw,DFE3Gamma,DFE3meanS_S99,DFE3beta)



layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE)) # four panel figure with one empty corner
#DFE1 plot first 
#locations 1, 4, 9 in list 
DFE1_matrix <- t(matrix(c(c(DFE1_truth[[1]][1],DFE1_truth[[2]][1],DFE1_truth[[3]][1],DFE1_truth[[4]][1]),
    c(discrete_classes[[1]][1],discrete_classes[[2]][1],discrete_classes[[3]][1],discrete_classes[[4]][1]),
    c(discrete_classes[[1]][4],discrete_classes[[2]][4],discrete_classes[[3]][4],discrete_classes[[4]][4]), 
    c(discrete_classes[[1]][7],discrete_classes[[2]][7],discrete_classes[[3]][7],discrete_classes[[4]][7])), 
    nr=4))
colnames(DFE1_matrix) <- c("f0", "f1", "f2", "f3")

barplot(DFE1_matrix, beside=T, ylab = "proportion mutations", main = "DFE1",
    col=c("black", "red","blue","green"))
legend("topleft", c("truth","selfing0","selfing50","selfing99"), pch=15, 
       col=c("black", "red","blue","green"), 
       bty="n")

# DFE2 
DFE2_matrix <- t(matrix(c(c(DFE2_truth[[1]][1],DFE2_truth[[2]][1],DFE2_truth[[3]][1],DFE2_truth[[4]][1]),
    c(discrete_classes[[1]][2],discrete_classes[[2]][2],discrete_classes[[3]][2],discrete_classes[[4]][2]),
    c(discrete_classes[[1]][5],discrete_classes[[2]][5],discrete_classes[[3]][5],discrete_classes[[4]][5]), 
    c(discrete_classes[[1]][8],discrete_classes[[2]][8],discrete_classes[[3]][8],discrete_classes[[4]][8])), 
    nr=4))
colnames(DFE2_matrix) <- c("f0", "f1", "f2", "f3")

barplot(DFE2_matrix, beside=T, ylab = "proportion mutations", main = "DFE2",
    col=c("black", "red","blue","green"))
legend("topleft", c("truth","selfing0","selfing50","selfing99"), pch=15, 
       col=c("black", "red","blue","green"), 
       bty="n")

# DFE3 
DFE3_matrix <- t(matrix(c(c(DFE3_truth[[1]][1],DFE3_truth[[2]][1],DFE3_truth[[3]][1],DFE3_truth[[4]][1]),
    c(discrete_classes[[1]][3],discrete_classes[[2]][3],discrete_classes[[3]][3],discrete_classes[[4]][3]),
    c(discrete_classes[[1]][6],discrete_classes[[2]][6],discrete_classes[[3]][6],discrete_classes[[4]][6]), 
    c(discrete_classes[[1]][9],discrete_classes[[2]][9],discrete_classes[[3]][9],discrete_classes[[4]][9])), 
    nr=4))
colnames(DFE3_matrix) <- c("f0", "f1", "f2", "f3")

barplot(DFE3_matrix, beside=T, ylab = "proportion mutations", main = "DFE3",
    col=c("black", "red","blue","green"))
legend("topleft", c("truth","selfing0","selfing50","selfing99"), pch=15, 
       col=c("black", "red","blue","green"), 
       bty="n")

x <- seq(0, 150, by=1)
#truth value gammas just to see what they look like
#plot(dgamma(x, 0.9, rate = 0.9/(5/200)), type = "l")
#abline(v = c(1,10,100), col = "red")
#plot(dgamma(x, 0.5, rate = 0.5/(50/200)), type = "l")
#abline(v = c(1,10,100), col = "red")
#plot(dgamma(x, 0.3, rate = 0.3/(1000/200)), type = "l")
#abline(v = c(1,10,100), col = "red")

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE)) # four panel figure with one empty corner
#DFE1, with new truth adjustments added
DFE1_matrix_mod <- t(matrix(c(c(DFE1_truth[[1]][1],DFE1_truth[[2]][1],DFE1_truth[[3]][1],DFE1_truth[[4]][1]),
    c(DFE1_truth_S0[[1]][1],DFE1_truth_S0[[2]][1],DFE1_truth_S0[[3]][1],DFE1_truth_S0[[4]][1]),
    c(discrete_classes[[1]][1],discrete_classes[[2]][1],discrete_classes[[3]][1],discrete_classes[[4]][1]),
    c(DFE1_truth_S50[[1]][1],DFE1_truth_S50[[2]][1],DFE1_truth_S50[[3]][1],DFE1_truth_S50[[4]][1]),
    c(discrete_classes[[1]][4],discrete_classes[[2]][4],discrete_classes[[3]][4],discrete_classes[[4]][4]),
    c(DFE1_truth_S99[[1]][1],DFE1_truth_S99[[2]][1],DFE1_truth_S99[[3]][1],DFE1_truth_S99[[4]][1]),
    c(discrete_classes[[1]][7],discrete_classes[[2]][7],discrete_classes[[3]][7],discrete_classes[[4]][7])), 
    nr=4))
colnames(DFE1_matrix_mod) <- c("f0", "f1", "f2", "f3")

barplot(DFE1_matrix_mod, beside=T, ylab = "proportion mutations", main = "DFE1",
    col=c("black", "#be5757","red","#41419b", "blue","#85be85", "green"))
legend("topleft", c("truth","selfing0 truth", "selfing0", 
    "selfing50 truth", "selfing50","selfing99 truth", "selfing99"), pch=15, 
       col=c("black", "#be5757","red","#41419b", "blue","#85be85", "green"), 
       bty="n")

#DFE2, with new truth adjustments added
DFE2_matrix_mod <- t(matrix(c(c(DFE2_truth[[1]][1],DFE2_truth[[2]][1],DFE2_truth[[3]][1],DFE2_truth[[4]][1]),
    c(DFE2_truth_S0[[1]][1],DFE2_truth_S0[[2]][1],DFE2_truth_S0[[3]][1],DFE2_truth_S0[[4]][1]),
    c(discrete_classes[[1]][2],discrete_classes[[2]][2],discrete_classes[[3]][2],discrete_classes[[4]][2]),
    c(DFE2_truth_S50[[1]][1],DFE2_truth_S50[[2]][1],DFE2_truth_S50[[3]][1],DFE2_truth_S50[[4]][1]),
    c(discrete_classes[[1]][5],discrete_classes[[2]][5],discrete_classes[[3]][5],discrete_classes[[4]][5]),
    c(DFE2_truth_S99[[1]][1],DFE2_truth_S99[[2]][1],DFE2_truth_S99[[3]][1],DFE2_truth_S99[[4]][1]),
    c(discrete_classes[[1]][8],discrete_classes[[2]][8],discrete_classes[[3]][8],discrete_classes[[4]][8])), 
    nr=4))
colnames(DFE2_matrix_mod) <- c("f0", "f1", "f2", "f3")

barplot(DFE2_matrix_mod, beside=T, ylab = "proportion mutations", main = "DFE2",
    col=c("black", "#be5757","red","#41419b", "blue","#85be85", "green"))
legend("topleft", c("truth","selfing0 truth", "selfing0", 
    "selfing50 truth", "selfing50","selfing99 truth", "selfing99"), pch=15, 
       col=c("black", "#be5757","red","#41419b", "blue","#85be85", "green"), 
       bty="n")

#DFE3, with new truth adjustments added
DFE3_matrix_mod <- t(matrix(c(c(DFE3_truth[[1]][1],DFE3_truth[[2]][1],DFE3_truth[[3]][1],DFE3_truth[[4]][1]),
    c(DFE3_truth_S0[[1]][1],DFE3_truth_S0[[2]][1],DFE3_truth_S0[[3]][1],DFE3_truth_S0[[4]][1]),
    c(discrete_classes[[1]][3],discrete_classes[[2]][3],discrete_classes[[3]][3],discrete_classes[[4]][3]),
    c(DFE3_truth_S50[[1]][1],DFE3_truth_S50[[2]][1],DFE3_truth_S50[[3]][1],DFE3_truth_S50[[4]][1]),
    c(discrete_classes[[1]][5],discrete_classes[[2]][5],discrete_classes[[3]][5],discrete_classes[[4]][5]),
    c(DFE3_truth_S99[[1]][1],DFE3_truth_S99[[2]][1],DFE3_truth_S99[[3]][1],DFE3_truth_S99[[4]][1]),
    c(discrete_classes[[1]][9],discrete_classes[[2]][9],discrete_classes[[3]][9],discrete_classes[[4]][9])), 
    nr=4))
colnames(DFE3_matrix_mod) <- c("f0", "f1", "f2", "f3")

barplot(DFE3_matrix_mod, beside=T, ylab = "proportion mutations", main = "DFE3",
    col=c("black", "#be5757","red","#41419b", "blue","#85be85", "green"))
legend("topleft", c("truth","selfing0 truth", "selfing0", 
    "selfing50 truth", "selfing50","selfing99 truth", "selfing99"), pch=15, 
       col=c("black", "#be5757","red","#41419b", "blue","#85be85", "green"), 
       bty="n")
