# now i will calculate B using Nordberg's equations for 10 neutral sites in my genome 
library(ggpubr)

file_content <- readLines("/nas/longleaf/home/adaigle/simulate_DFEs/genome.txt")

g3_elements <- list()

# Process each line of the file
for (line in file_content) {
  # Extract the genomic element type, start position, and end position from each line
  line_elements <- strsplit(line, ",")[[1]]
  element_type <- gsub("initializeGenomicElement\\(|\\s+", "", line_elements[1])
  # Check if the element is of type g3
  if (element_type == "g3") {
    start_position <- as.numeric(gsub("\\s+", "", line_elements[2]))
    end_position <- as.numeric(gsub("\\s+|\\);", "", line_elements[3]))
    
    # Create a sequence of numbers from start to end
    element_numbers <- seq(start_position, end_position)
    
    # Append the element numbers to the list
    g3_elements[[length(g3_elements) + 1]] <- element_numbers
  }
}

computational_B <- function(meanGamma,beta,pself) {
# Create the nordborg table for calculations
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- ifelse(nordborg_calcs$Assignment == "selected", 
  rgamma(nrow(nordborg_calcs), shape = s_shape, rate = s_rate), 0)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)

Bs <- numeric()
neutral_sites_list <- numeric()
for (x in seq(1,1000, by = 1)) {
neutral_position  <- neutral_sites[sample(1:nrow(neutral_sites), 1), ]

neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 *100#mut rate
rec_rate <- 3.12*1e-8 *100 # rec rate
 
selected_table$ti <- selected_table$s * 0.5
selected_table$haldanesr <- 1- exp((-2*rec_rate*selected_table$Distance)/2)

#selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (selected_table$haldanesr*(1-F))*(1-selected_table$ti))^2
selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2
#selected_table$nordborg <- (ui*selected_table$ti) / 2*(selected_table$ti + (rec_rate*selected_table$Distance*(1-F)))^2

B <- exp(-sum(selected_table$nordborg))
  # Store the results
  Bs <- c(Bs, B)
}
genome_B <- data.frame(Bs, neutral_sites_list, distance)
return(genome_B)
}

tableDFE1 <- computational_B(5,0.9,0)
mean(tableDFE1$Bs)
sd(tableDFE1$Bs)
plot(tableDFE1$neutral_sites_list, tableDFE1$Bs, ylim=c(0,1))

tableDFE2 <- computational_B(50,0.5,0)
mean(tableDFE2$Bs)
sd(tableDFE2$Bs)
plot(tableDFE2$neutral_sites_list, tableDFE2$Bs, ylim=c(0.8,1))

tableDFE3 <- computational_B(1000,0.3,0)
mean(tableDFE3$Bs)
sd(tableDFE3$Bs)
plot(tableDFE3$neutral_sites_list, tableDFE2$Bs, ylim=c(0.8,1))

#averages <- numeric()
#standard_deviations <- numeric()
#for (x in seq(0.01,4000, by = 100)) {
#  DFE2_results <- lapply(1:10, function(y) computational_B(x, 0.5,0))
#  DFE2_results <- unlist(DFE2_results)
#  
#  # Calculate the average and standard deviation
#  average <- median(DFE2_results)
#  sd_value <- sd(DFE2_results)
#  
#  # Store the results
#  averages <- c(averages, average)
#  standard_deviations <- c(standard_deviations, sd_value)
#}
#
## Plot the averages and standard deviations
#plot(seq(0.01, 4000, by = 100), averages, type = "l", xlab = "Mean Gamma", ylab = "average B")
#lines(seq(0.01, 4000, by = 100), standard_deviations, col = "red")
#legend("topright", legend = c("Average", "Standard Deviation"), col = c("black", "red"), lty = 1)



computational_B_center_end_comparison <- function(meanGamma,beta,pself) {
# Create the nordborg table for calculations
#using first 10kb, middle 10kb and last 10kb
#since neutral exonic sites are .0375 prop of sites, this is 375 sites for each region
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- ifelse(nordborg_calcs$Assignment == "selected", 
  rgamma(nrow(nordborg_calcs), shape = s_shape, rate = s_rate), 0)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)

Bs <- numeric()
neutral_sites_list <- numeric()
len <- nrow(neutral_sites)
sites <- c(seq(1,375, by = 1), seq((len/2-187),(len/2+188), by = 1), seq(len-375,len, by = 1))
for (x in sites) {
neutral_position  <- neutral_sites[x, ]

neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 * 100#mut rate
rec_rate <- 3.12*1e-8 * 100 # rec rate

 
selected_table$ti <- selected_table$s * 0.5

selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2

B <- exp(-sum(selected_table$nordborg))
  # Store the results
  Bs <- c(Bs, B)
}
genome_B <- data.frame(Bs, neutral_sites_list)
return(genome_B)
}

#tableDFE2 <- computational_B_center_end_comparison(50,0.5,.00001)
#plot(tableDFE2$neutral_sites_list, tableDFE2$Bs, ylim=c(0.8,1))
#
#tableDFE3 <- computational_B_center_end_comparison(1000,0.3, .99)
#plot(tableDFE3$neutral_sites_list, tableDFE3$Bs, ylim=c(0,1))



tableDFE20 <- computational_B(50,0.5,0)

tableDFE20001 <- computational_B(50,0.5,.0001)

tableDFE250 <- computational_B(50,0.5,.50)

tableDFE280 <- computational_B(50,0.5,.80)
tableDFE299 <- computational_B(50,0.5,.90)

tableDFE299 <- computational_B(50,0.5,.99)

p1 <- plot(tableDFE20$neutral_sites_list, tableDFE20$Bs, ylim=c(0,1))
p2 <- plot(tableDFE20001$neutral_sites_list, tableDFE20001$Bs, ylim=c(0,1))
p3 <- plot(tableDFE250$neutral_sites_list, tableDFE250$Bs, ylim=c(0,1))
p4 <- plot(tableDFE280$neutral_sites_list, tableDFE280$Bs, ylim=c(0,1))
p5 <- plot(tableDFE299$neutral_sites_list, tableDFE299$Bs, ylim=c(0,1))

ggarrange(p1, p2, p3, p4, p5, 
                    labels = c("0", "0.0001", "0.5", "0.8", "0.99"),
                    font.label = list(size = 24, color = "black", face = "bold", family = NULL),
                    ncol = 2, nrow = 3,
                    common.legend = TRUE, legend = "right", vjust=1)

tableDFE3_0 <- computational_B(1000,0.3,0)
plot(tableDFE3_0$neutral_sites_list, tableDFE3_0$Bs, ylim=c(0,1))
mean(tableDFE3_0$Bs)
mean(tableDFE20$Bs)

ui <- 3.3e-9 #mut rate
rec_rate <- 3.12*1e-8 # rec rate

get_s_dist <- function(meanGamma,beta,pself) {
# For one neutral site get B by distance
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- ifelse(nordborg_calcs$Assignment == "selected", 
  rgamma(nrow(nordborg_calcs), shape = s_shape, rate = s_rate), 0)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)


neutral_position  <- neutral_sites[sample(1:nrow(neutral_sites), 1), ]

#neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 *100#mut rate
rec_rate <- 3.12*1e-8 *100 # rec rate
 
selected_table$ti <- selected_table$s * 0.5
selected_table$haldanesr <- 1- exp((-2*rec_rate*selected_table$Distance)/2)

#selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (selected_table$haldanesr*(1-F))*(1-selected_table$ti))^2
selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2
#selected_table$nordborg <- (ui*selected_table$ti) / 2*(selected_table$ti + (rec_rate*selected_table$Distance*(1-F)))^2

selected_table$B <- exp(-(selected_table$nordborg))
selected_table <- selected_table[selected_table$Distance <= 20000, ]
selected_table_subset <- selected_table[sample(1:nrow(selected_table), 1000), ]
return(selected_table_subset)
}

s_dist_DFE2 <- get_s_dist(50,0.5,0)
plot(s_dist_DFE2$Distance, s_dist_DFE2$B, ylim=c(0.994,1))

s_dist_DFE2_99 <- get_s_dist(50,0.5,.99)
plot(s_dist_DFE2_99$Distance, s_dist_DFE2_99$B, ylim=c(0.994,1))






computational_B_start <- function(meanGamma,beta,pself) {
# Create the nordborg table for calculations
#using first 10kb, middle 10kb and last 10kb
#since neutral exonic sites are .0375 prop of sites, this is 375 sites for each region
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- ifelse(nordborg_calcs$Assignment == "selected", 
  rgamma(nrow(nordborg_calcs), shape = s_shape, rate = s_rate), 0)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)

Bs <- numeric()
neutral_sites_list <- numeric()
len <- nrow(neutral_sites)
sites <- c(seq(1,375, by = 1))
for (x in sites) {
neutral_position  <- neutral_sites[x, ]

neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 * 100#mut rate
rec_rate <- 3.12*1e-8 * 100 # rec rate

 
selected_table$ti <- selected_table$s * 0.5

selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2

B <- exp(-sum(selected_table$nordborg))
  # Store the results
  Bs <- c(Bs, B)
}
genome_B <- data.frame(Bs, neutral_sites_list)
return(genome_B)
}


computational_B_middle <- function(meanGamma,beta,pself) {
# Create the nordborg table for calculations
#using first 10kb, middle 10kb and last 10kb
#since neutral exonic sites are .0375 prop of sites, this is 375 sites for each region
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- ifelse(nordborg_calcs$Assignment == "selected", 
  rgamma(nrow(nordborg_calcs), shape = s_shape, rate = s_rate), 0)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)

Bs <- numeric()
neutral_sites_list <- numeric()
len <- nrow(neutral_sites)
sites <- c(seq((len/2-187),(len/2+188), by = 1))
for (x in sites) {
neutral_position  <- neutral_sites[x, ]

neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 * 100#mut rate
rec_rate <- 3.12*1e-8 * 100 # rec rate

 
selected_table$ti <- selected_table$s * 0.5

selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2

B <- exp(-sum(selected_table$nordborg))
  # Store the results
  Bs <- c(Bs, B)
}
genome_B <- data.frame(Bs, neutral_sites_list)
return(genome_B)
}


computational_B_end <- function(meanGamma,beta,pself) {
# Create the nordborg table for calculations
#using first 10kb, middle 10kb and last 10kb
#since neutral exonic sites are .0375 prop of sites, this is 375 sites for each region
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- ifelse(nordborg_calcs$Assignment == "selected", 
  rgamma(nrow(nordborg_calcs), shape = s_shape, rate = s_rate), 0)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)

Bs <- numeric()
neutral_sites_list <- numeric()
len <- nrow(neutral_sites)
sites <- c(seq(1,375, by = 1), seq((len/2-187),(len/2+188), by = 1), seq(len-375,len, by = 1))
for (x in sites) {
neutral_position  <- neutral_sites[x, ]

neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 * 100#mut rate
rec_rate <- 3.12*1e-8 * 100 # rec rate

 
selected_table$ti <- selected_table$s * 0.5

selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2

B <- exp(-sum(selected_table$nordborg))
  # Store the results
  Bs <- c(Bs, B)
}
genome_B <- data.frame(Bs, neutral_sites_list)
return(genome_B)
}

tableDFE2_start <- computational_B_start(5,0.9,0)


tableDFE2_middle <- computational_B_middle(5,0.9,0)


tableDFE2_end <- computational_B_end(5,0.9,0)
mean(tableDFE2_start$Bs)
sd(tableDFE2_start$Bs)
mean(tableDFE2_middle$Bs)
sd(tableDFE2_middle$Bs)
mean(tableDFE2_end$Bs)
sd(tableDFE2_end$Bs)

tableDFE3_start <- computational_B_start(1000,0.3,0)
tableDFE3_middle <- computational_B_middle(1000,0.3,0)
tableDFE3_end <- computational_B_end(1000,0.3,0)

mean(tableDFE3_start$Bs)
sd(tableDFE3_start$Bs)
mean(tableDFE3_middle$Bs)
sd(tableDFE3_middle$Bs)
mean(tableDFE3_end$Bs)
sd(tableDFE3_end$Bs)

tableDFE20 <- computational_B(50,0.5,0)
tableDFE20 <- computational_B(1000,0.3,0)


computational_B_constant <- function(meanGamma,beta,pself) {
# Create the nordborg table for calculations
F <- pself / (2 - pself)

nordborg_calcs <- data.frame(Number = unlist(g3_elements))

nordborg_calcs$Assignment <- sample(c("selected", "neutral"), size = nrow(nordborg_calcs), replace = TRUE)
meanS <- meanGamma/(2*5000)
s_shape <- beta
s_rate <- s_shape/abs(meanS)

nordborg_calcs$s <- meanGamma / (2*5000)

neutral_sites <- subset(nordborg_calcs, Assignment == "neutral")
#center_index <- ceiling(nrow(neutral_sites) / 2)

Bs <- numeric()
neutral_sites_list <- numeric()
for (x in seq(1,1000, by = 1)) {
neutral_position  <- neutral_sites[sample(1:nrow(neutral_sites), 1), ]

neutral_sites_list <- c(neutral_sites_list, neutral_position$Number)
selected_table <- subset(nordborg_calcs, Assignment == "selected")

# Calculate the distance of each selected site from the neutral site
selected_table$Distance <- abs(selected_table$Number - neutral_position$Number)

ui <- 3.3e-9 *100#mut rate
rec_rate <- 3.12*1e-8 *100 # rec rate
 
selected_table$ti <- selected_table$s * 0.5
selected_table$haldanesr <- 1- exp((-2*rec_rate*selected_table$Distance)/2)

#selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (selected_table$haldanesr*(1-F))*(1-selected_table$ti))^2
selected_table$nordborg <- (ui*selected_table$ti) / (selected_table$ti + (rec_rate*selected_table$Distance*(1-F))*(1-selected_table$ti))^2
#selected_table$nordborg <- (ui*selected_table$ti) / 2*(selected_table$ti + (rec_rate*selected_table$Distance*(1-F)))^2

B <- exp(-sum(selected_table$nordborg))
  # Store the results
  Bs <- c(Bs, B)
}
genome_B <- data.frame(Bs, neutral_sites_list, distance)
return(genome_B)
}

tableDFE3_const <- computational_B_constant(1000,0.3,0)
mean(tableDFE3_const$Bs)
