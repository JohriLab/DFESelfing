rm(list=ls())
library(RColorBrewer)
library(tidyverse)
source("/nas/longleaf/home/adaigle/DFESelfing/calculate_pi.r")
library(reshape2)

dfealpha_dir <- "/nas/longleaf/home/adaigle/DFESelfing/"
dfealpha_output_dirs <- paste(paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_output"), sep = ""),
    "/", dir(
    paste(dfealpha_dir, dir(dfealpha_dir, pattern = "DFE_alpha_output"), sep = ""), pattern = "selected")
    , sep = "")


#table of results with columns signifying selfing%, DFE, and output
dfealpha_colnames <- as.data.frame(c("N1", "N2", "t2", "Nw", "b", "Es", "f0","L"))
dfealpha_raw_results <- tibble(
    name = list.files(path = dfealpha_output_dirs, pattern = "est_dfe.out"),
    fullpath = paste(dfealpha_output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=DFE_alpha_output_)\\d+") %>% if_else(is.na(.), "0", .), #0's were empty
    matchname = str_extract(fullpath, "(DFE)\\d+(output)\\d"), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    data = lapply(fullpath, function(x) 
        as.data.frame(t(data.frame(matrix(unlist(strsplit(readLines(x), " ")), ncol = 2, byrow = TRUE)))) %>%
            rownames_to_column() %>%
            `colnames<-`(.[1,]) %>%
            .[-1,] %>%
            `rownames<-`(NULL))
)

df <- dfealpha_raw_results$data %>% bind_rows()
df <- df %>% mutate(row_number = row_number())
dfealpha_raw_results <- dfealpha_raw_results  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
dfealpha_raw_results <- dfealpha_raw_results %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame %>% 
    mutate(across(.fns = ~ if(all(!is.na(as.numeric(.x)))) as.numeric(.x) else .x)) %>%
    mutate(gamma = -2*Nw*Es) %>%
    rename(f0_fromoutput = f0) #f0 is a dfealpha output, but I also need to name a class f0 later

DFE_proportions_dfe_alpha <- function(meanGamma,beta, Nw) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    meanS <- meanGamma/(2*Nw)

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}
true_gammas <- list(
  DFE1 = 5,
  DFE2 = 50,
  DFE3 = 1000
)
true_betas <- list(
  DFE1 = 0.9,
  DFE2 = 0.5,
  DFE3 = 0.3
)

dfealpha_raw_results_wtruth <- dfealpha_raw_results %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

#now i run the class generator function on all outputs
dfealpha_raw_results_wclasses <- dfealpha_raw_results_wtruth %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) %>%
  rename(z0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(z1 = f1) %>% 
  rename(z2 = f2) %>% 
  rename(z3 = f3) %>% 
  #rename( c(t0,t1,t2,t3) = c(f0,f1,f2,f3)) %>%
  mutate(output = pmap(list(gamma, b, Nw), DFE_proportions_dfe_alpha)) %>%
  unnest_wider(output) %>%
  mutate(across(.fns = ~ if(all(!is.na(as.numeric(.x)))) as.numeric(.x) else .x))

dfealpha_summary <- dfealpha_raw_results_wclasses %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd))) %>%
    rename(z0 = z0_avg) %>%  #had to rename columns bc i rerun the function
    rename(z1 = z1_avg) %>% 
    rename(z2 = z2_avg) %>% 
    rename(z3 = z3_avg) %>%
    rename(f0 = f0_avg) %>%  #had to rename columns bc i rerun the function
    rename(f1 = f1_avg) %>% 
    rename(f2 = f2_avg) %>% 
    rename(f3 = f3_avg) 

  # convert the data frame to a tidy format
dfealpha_tidy <- gather(dfealpha_summary, 
    key = "generation", value = "value", c(z0,z1,z2,z3)) %>%
    mutate(generation = recode(generation,
    z0 = 'f0', 
    z1 = 'f1', 
    z2 = 'f2', 
    z3 = 'f3')) %>%
    mutate(selfing = as.character(selfing))

##plotting done at bottom bc I need grapes 

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
selfing80 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output_80/")
selfing99 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output_100/")

selfing90 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output_90/")
selfing95 <- ("/nas/longleaf/home/adaigle/DFESelfing/DFE_alpha_output_95/")



#run function to create dataframes for each selfing percentage
means0 <- find_means(selfing0, DFE_list, replicates)
means50 <- find_means(selfing50, DFE_list, replicates)
means80 <- find_means(selfing80, DFE_list, replicates)
means99 <- find_means(selfing99, DFE_list, replicates)

means90 <- find_means(selfing90, DFE_list, replicates)
means95 <- find_means(selfing95, DFE_list, replicates)


#betas and gammas used to generate DFE1,2,3
truebeta <- c(0.9,0.5,0.3)
truegamma <- c(5,50,1000)

#makes a summary table showing accuracy of beta and gamma predictions
gammabetatable <- data.frame(truebeta, means0$b, means0$betasd, means50$b, means50$betasd,
means80$b, means80$betasd, means90$b, means90$betasd, means95$b, means95$betasd, 
means99$b, means99$betasd, truegamma, -as.numeric(means0$gamma), means0$gammasd, 
-as.numeric(means50$gamma), means50$gammasd, -as.numeric(means80$gamma), means80$gammasd,
-as.numeric(means90$gamma), means90$gammasd,-as.numeric(means95$gamma), means95$gammasd,
-as.numeric(means99$gamma), means99$gammasd)

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

Nw80 <- as.numeric(means80$Nw)
#meanS <- as.numeric(args[2]) # shape = mean * rate
meanGamma80 <- -as.numeric(means80$gamma)
meanS80 <- meanGamma80/(2.0*Nw80)
beta80 <- as.numeric(means80$b) #shape parameter

#using 2NwS for gamma for now, wasn't quite sure what to do. 
#vector of lists, f0, f1,f2, f4 (discrete mutation classes)
#within each class is selfing levels(DFE1-3)
#so 0selfing DFE1,2,3 ; 50 selfing DFE1-3 etc.
dfe80<-DFE_proportions(Nw80,meanGamma80,meanS80,beta80)

#---------------quick dfe90 and 95 stuff
Nw90 <- as.numeric(means90$Nw)
#meanS <- as.numeric(args[2]) # shape = mean * rate
meanGamma90 <- -as.numeric(means90$gamma)
meanS90 <- meanGamma90/(2.0*Nw90)
beta90 <- as.numeric(means90$b) #shape parameter

#using 2NwS for gamma for now, wasn't quite sure what to do. 
#vector of lists, f0, f1,f2, f4 (discrete mutation classes)
#within each class is selfing levels(DFE1-3)
#so 0selfing DFE1,2,3 ; 50 selfing DFE1-3 etc.
dfe90<-DFE_proportions(Nw90,meanGamma90,meanS90,beta90)

Nw95 <- as.numeric(means95$Nw)
#meanS <- as.numeric(args[2]) # shape = mean * rate
meanGamma95 <- -as.numeric(means95$gamma)
meanS95 <- meanGamma95/(2.0*Nw95)
beta95 <- as.numeric(means95$b) #shape parameter

#using 2NwS for gamma for now, wasn't quite sure what to do. 
#vector of lists, f0, f1,f2, f4 (discrete mutation classes)
#within each class is selfing levels(DFE1-3)
#so 0selfing DFE1,2,3 ; 50 selfing DFE1-3 etc.
dfe95<-DFE_proportions(Nw95,meanGamma95,meanS95,beta95)


layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE)) # four panel figure with one empty corner
#DFE1 plot first 
#locations 1, 4, 9 in list 
DFE1_matrix <- t(matrix(c(c(DFE1_truth[[1]][1],DFE1_truth[[2]][1],DFE1_truth[[3]][1],DFE1_truth[[4]][1]),
    c(discrete_classes[[1]][1],discrete_classes[[2]][1],discrete_classes[[3]][1],discrete_classes[[4]][1]),
    c(discrete_classes[[1]][4],discrete_classes[[2]][4],discrete_classes[[3]][4],discrete_classes[[4]][4]),
    c(dfe80[[1]][1],dfe80[[2]][1],dfe80[[3]][1],dfe80[[4]][1]),
    c(dfe90[[1]][1],dfe90[[2]][1],dfe90[[3]][1],dfe90[[4]][1]),
    c(dfe95[[1]][1],dfe95[[2]][1],dfe95[[3]][1],dfe95[[4]][1]),
    c(discrete_classes[[1]][7],discrete_classes[[2]][7],discrete_classes[[3]][7],discrete_classes[[4]][7])), 
    nr=4))
colnames(DFE1_matrix) <- c("f0", "f1", "f2", "f3")

barplot(DFE1_matrix, beside=T, ylab = "proportion mutations", main = "DFE1",
    col=c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"))
legend("topleft", c("truth","selfing0","selfing50","selfing80", "selfing90", "selfing95", "selfing99"), pch=15, 
       col=c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"), 
       bty="n")

# DFE2 
DFE2_matrix <- t(matrix(c(c(DFE2_truth[[1]][1],DFE2_truth[[2]][1],DFE2_truth[[3]][1],DFE2_truth[[4]][1]),
    c(discrete_classes[[1]][2],discrete_classes[[2]][2],discrete_classes[[3]][2],discrete_classes[[4]][2]),
    c(discrete_classes[[1]][5],discrete_classes[[2]][5],discrete_classes[[3]][5],discrete_classes[[4]][5]), 
    c(dfe80[[1]][2],dfe80[[2]][2],dfe80[[3]][2],dfe80[[4]][2]), 
    c(dfe90[[1]][2],dfe90[[2]][2],dfe90[[3]][2],dfe90[[4]][2]),
    c(dfe95[[1]][2],dfe95[[2]][2],dfe95[[3]][2],dfe95[[4]][2]),
    c(discrete_classes[[1]][8],discrete_classes[[2]][8],discrete_classes[[3]][8],discrete_classes[[4]][8])), 
    nr=4))
colnames(DFE2_matrix) <- c("f0", "f1", "f2", "f3")

barplot(DFE2_matrix, beside=T, ylab = "proportion mutations", main = "DFE2",
    col=c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"))
legend("topleft", c("truth","selfing0","selfing50","selfing80", "selfing90", "selfing95", "selfing99"), pch=15, 
       col=c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"), 
       bty="n")

# DFE3 
DFE3_matrix <- t(matrix(c(c(DFE3_truth[[1]][1],DFE3_truth[[2]][1],DFE3_truth[[3]][1],DFE3_truth[[4]][1]),
    c(discrete_classes[[1]][3],discrete_classes[[2]][3],discrete_classes[[3]][3],discrete_classes[[4]][3]),
    c(discrete_classes[[1]][6],discrete_classes[[2]][6],discrete_classes[[3]][6],discrete_classes[[4]][6]),
    c(dfe80[[1]][3],dfe80[[2]][3],dfe80[[3]][3],dfe80[[4]][3]),
    c(dfe90[[1]][3],dfe90[[2]][3],dfe90[[3]][3],dfe90[[4]][3]),
    c(dfe95[[1]][3],dfe95[[2]][3],dfe95[[3]][3],dfe95[[4]][3]),
    c(discrete_classes[[1]][9],discrete_classes[[2]][9],discrete_classes[[3]][9],discrete_classes[[4]][9])), 
    nr=4))
colnames(DFE3_matrix) <- c("f0", "f1", "f2", "f3")

barplot(DFE3_matrix, beside=T, ylab = "proportion mutations", main = "DFE3",
    col=c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"))
legend("topleft", c("truth","selfing0","selfing50","selfing80", "selfing90", "selfing95", "selfing99"), pch=15, 
       col=c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"), 
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
    c(discrete_classes[[1]][6],discrete_classes[[2]][6],discrete_classes[[3]][6],discrete_classes[[4]][6]),
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


DFE3_matrix_mod <- t(matrix(c(c(DFE3_truth[[1]][1],DFE3_truth[[2]][1],DFE3_truth[[3]][1],DFE3_truth[[4]][1]),
    c(DFE3_truth_S0[[1]][1],DFE3_truth_S0[[2]][1],DFE3_truth_S0[[3]][1],DFE3_truth_S0[[4]][1]),
    c(discrete_classes[[1]][3],discrete_classes[[2]][3],discrete_classes[[3]][3],discrete_classes[[4]][3]),
    c(DFE3_truth_S50[[1]][1],DFE3_truth_S50[[2]][1],DFE3_truth_S50[[3]][1],DFE3_truth_S50[[4]][1]),
    c(discrete_classes[[1]][6],discrete_classes[[2]][6],discrete_classes[[3]][6],discrete_classes[[4]][6]),
    c(DFE3_truth_S99[[1]][1],DFE3_truth_S99[[2]][1],DFE3_truth_S99[[3]][1],DFE3_truth_S99[[4]][1]),
    c(discrete_classes[[1]][9],discrete_classes[[2]][9],discrete_classes[[3]][9],discrete_classes[[4]][9])), 
    nr=4))


# start grapes stuff, keeping separate for now
rm(list=ls())
library(tidyverse)

#read in all selfing %, dfes, and experiments
grapes_dir <- "/nas/longleaf/home/adaigle/DFESelfing/grapes/"
output_dirs <- paste(grapes_dir, dir(grapes_dir, pattern = "output"), "/all", sep = "")

#table of results with columns signifying selfing%, DFE, and output
grapes_raw_results <- tibble(
    name = list.files(path = output_dirs, pattern = ".csv$"),
    matchname = sub(".txt.csv","",name), #will be useful for comparing results later, or combining with other programs
    DFE = str_extract(matchname, "(DFE)\\d+"),
    fullpath = paste(output_dirs, "/", name,sep=""),
    selfing = str_extract(fullpath, "(?<=grapes_output_)\\d+"),
    data = lapply(fullpath,read.csv)
)

#now I want to clean this up because it has too much going on
# I will keep just the gammazero row in each DF
grapes_gammazero_raw_results <- grapes_raw_results %>% 
    mutate(data = map(data, ~ { #start of a mapping funciton, which is applied to each df in column
    df_filtered <- .x %>% filter(model == 'GammaZero') %>% #select only GammaZero row
    select(1:52) %>%# removing prediction columns for now
    select_if(~ all(!is.na(.))) # remove NA columns
    return(df_filtered)
}))

#now i summarize the outputs by averaging them and finding the standard deviation
#first I turn my table of dataframes into a dataframe, because my df's only had one row 
# Extract the column of data frames and bind them into a single data frame
df <- grapes_gammazero_raw_results$data %>% bind_rows()

# Add a column to the data frame that contains the row numbers
df <- df %>% mutate(row_number = row_number())
grapes_gammazero_raw_results <- grapes_gammazero_raw_results  %>% 
    mutate(row_number = row_number())

# Join the original table with the new data frame using the row numbers as the key
grapes_gammazero_raw_results <- grapes_gammazero_raw_results %>% 
    left_join(df, by = "row_number") %>% 
    select(-data, -row_number) %>%
    as.data.frame

# summarize using group_by and summarize
# I originally tried using table of dataframes, but got weird
# for now I'll take one experiment at a time

grapes_gammazero_summary <- grapes_gammazero_raw_results %>% 
    group_by(selfing, DFE) %>%
    summarize(across(where(is.numeric), list(avg = mean, sd = sd)))

true_gammas <- list(
  DFE1 = 5,
  DFE2 = 50,
  DFE3 = 1000
)
true_betas <- list(
  DFE1 = 0.9,
  DFE2 = 0.5,
  DFE3 = 0.3
)

grapes_gammazero_simple_summary <- 
    grapes_gammazero_summary[
        c('DFE','selfing', 'GammaZero.negGmean_avg', 
            'GammaZero.negGmean_sd', 'GammaZero.negGshape_avg',
            'GammaZero.negGshape_sd')] %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))



grapes_gammazero_simple_summary <- grapes_gammazero_simple_summary %>%
    group_by(DFE) %>%
    mutate(true_mean = unlist(true_gammas[DFE])) %>%
    mutate(true_shape = unlist(true_betas[DFE]))

#now make dfe class prediction for each gmean(gamma) and shape(b)
DFE_proportions_grapes <- function(meanGamma,beta) { 
    # modified function to calc vector of dfe classes given gamma and beta
    # assumes Nw is 100 bc grapes doesn't have two step
    # need to confirm it has no modified Ne
    Nw <- 100 # assuming this is what we should do? 
    meanS <- meanGamma/(2*Nw)

    # this code adapted from get_dfe_class_proportions
    s_shape <- beta
    s_rate <- s_shape/abs(meanS)

    #print("gamma distribution with shape parameter, alpha:")
    #print(s_shape)
    #print("gamma distribution with rate parameter, beta:")
    #print(s_rate)
    #these statements refer to the wikipedia way of describing gamma
    #https://en.wikipedia.org/wiki/Gamma_distribution
    #not how r categorizes them!

    #print ("DFE classes in terms of Nws: ")
    #x1 <- 1/(Nw)
    #x10 <- 10/(Nw)
    #x100 <- 100/(Nw)
    #f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    #f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    #f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    #f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))

    #print ("DFE classes in terms of 2Nws: ")
    x1 <- 1/(2*Nw)
    x10 <- 10/(2*Nw)
    x100 <- 100/(2*Nw)
    f0 <- pgamma(x1, shape=s_shape, rate=s_rate)
    f1 <- pgamma(x10, shape=s_shape, rate=s_rate) - f0
    f2 <- pgamma(x100, shape=s_shape, rate=s_rate) - f1 - f0
    f3 <- 1.0 - f2 - f1 - f0
    #print(paste(f0, f1, f2, f3, sep=", "))
    #return(c(list(f0), list(f1), list(f2), list(f3)))
    return(c(f0 = f0, f1 = f1, f2 = f2, f3 = f3))
}

# run DFE_proportions funciton and save as new columns in df 
grapes_gammazero_simple_summary <- grapes_gammazero_simple_summary %>%
  mutate(output = map2(GammaZero.negGmean_avg, GammaZero.negGshape_avg, DFE_proportions_grapes)) %>%
  unnest_wider(output) %>%
  rename(t0 = f0) %>%  #had to rename columns bc i rerun the function
  rename(t1 = f1) %>% 
  rename(t2 = f2) %>% 
  rename(t3 = f3) %>% 
  #rename( c(t0,t1,t2,t3) = c(f0,f1,f2,f3)) %>%
  mutate(output = map2(true_mean, true_shape, DFE_proportions_grapes)) %>%
  unnest_wider(output) 

# convert the data frame to a tidy format
df_tidy <- gather(grapes_gammazero_simple_summary, 
    key = "generation", value = "value", c(t0:t3)) %>%
    mutate(generation = recode(generation,
    t0 = 'f0', 
    t1 = 'f1', 
    t2 = 'f2', 
    t3 = 'f3')) 

df_true_tidy <- gather(grapes_gammazero_simple_summary, 
    key = "generation", value = "value", c(f0:f3)) %>%
    select(1,13:14) %>% distinct() %>%
    mutate(selfing = "True0")

# create a bar plot with separate panels for each DFE
ggplot(bind_rows(df_true_tidy,df_tidy) , aes(x = generation, y = value, fill = selfing)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") 


# Define the order of the selfing levels
selfing_order <- c("True0", 0, 50, 80, 90, 95, 99)

# Create the grouped bar chart with custom selfing order
ggplot(bind_rows(df_true_tidy,df_tidy), aes(x = generation, y = value, fill = factor(selfing, levels = c("True0", "0", "50", "80", "90", "95", "99")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  scale_fill_manual(values = c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"))

  
dfealpha_for_plotting <- bind_rows(df_true_tidy,dfealpha_tidy) %>%
    select(1:4,32,34,36,38) %>%
    melt() %>% 
    mutate(value = ifelse(is.na(value), 0, value))

# filter sd's
dfealpha_for_plotting <- dfealpha_for_plotting %>% 
  filter(variable == "value" | paste0(generation, "_sd") == variable)

dfealpha_for_plotting$variable <- ifelse(grepl("_sd", dfealpha_for_plotting$variable), "sd", "value")  

voodoo <- pivot_wider(dfealpha_for_plotting, id_cols = c("generation","DFE","selfing"), names_from = "variable", values_from = "value")
# Define the order of the selfing levels
selfing_order <- c("True0", 0, 50, 80, 90, 95, 100)

# Create the grouped bar chart with custom selfing order
ggplot(voodoo, aes(x = generation, y = value, fill = factor(selfing, levels = c("True0", "0", "50", "80", "90", "95", "100")))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ DFE, nrow = 2) +
  labs(title = "DFEalpha", x = "Mutation Class (least to most deleterious)", y = "proportion of mutations", fill = "Selfing %") +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("black", "red","blue","purple", "goldenrod4", "hotpink", "green"))
#now I take the average, getting avg and sd for my classes
#I could also probably just ggplot but I think we will always 
#be using the averages for now

# create example data
df <- data.frame(generation = rep(c("f0", "f1", "f2", "f3"), each = 5),
                 variable = rep(c("value", "f0_sd", "f1_sd", "f2_sd", "f3_sd"), times = 4))



  paste0()