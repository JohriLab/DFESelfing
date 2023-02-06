#This is to get the proportion of f0, f1, f2 and f3 in case of diploid and haploid:

args = commandArgs(trailingOnly=TRUE)
Nw <- as.numeric(args[1])
#meanS <- as.numeric(args[2]) # shape = mean * rate
meanGamma <- as.numeric(args[2])
meanS <- meanGamma/(2.0*Nw)
beta <- as.numeric(args[3]) #shape parameter

s_shape <- beta
s_rate <- s_shape/abs(meanS)

print("gamma distribution with shape parameter, alpha:")
print(s_shape)
print("gamma distribution with rate parameter, beta:")
print(s_rate)

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

print("done")

