find_subpopulation_size <- function(migration_rate, demes) {
    # find the correct subpop size in a fully-connected island model 
    # of pop structure for a given number of demes 
    # here the goal is to make the Nem (Ne of the total metapopulation) 
    #equal to 5000
    # eq 7.3c on pg 318 of Charlesworth and Charlesworth
    subpop_Ne <- ( (5000 - ((demes-1)^2) / 4*demes*migration_rate) ) / demes

    return(subpop_Ne)
}

#what migration rates should i use???
find_migration_rate <- function(subpop_Ne, Nem) {
    migration_rate <- Nem/subpop_Ne
    return(migration_rate)
}

find_subpopulation_size2 <- function(migration_rate, demes) {
    # find the correct subpop size in a fully-connected island model 
    # of pop structure for a given number of demes 
    # here the goal is to make the Nem (Ne of the total metapopulation) 
    #equal to 5000
    # eq 7.3c on pg 318 of Charlesworth and Charlesworth
    subpop_Ne <- (-5000+((demes-1)^2 / (4*demes*migration_rate))) / (-demes )

    return(subpop_Ne)
}


find_subpopulation_size3 <- function(migration_rate, demes) {
    # find the correct subpop size in a fully-connected island model 
    # of pop structure for a given number of demes 
    # here the goal is to make the Nem (Ne of the total metapopulation) 
    #equal to 5000
    # eq 7.3c on pg 318 of Charlesworth and Charlesworth
    subpop_Ne <- (20000*demes*migration_rate - (demes-1)^2 ) / (4*migration_rate*demes^2)

    return(subpop_Ne)
}

fst <- function(Ne, m, d) {
    return(1/(1+((4*Ne*m*d^2)/(d-1)^2)))
}
fst(385, 0.00026, 5)
fst(5000, 0.00026, 5)
fst(385*5, 0.00026, 5)

fst2 <- function(Ne, m, d) {
    return(1/(1+4*Ne*m))
}
fst2(385, 0.00026, 5)
fst2(385*5, 0.00026, 5)
fst2(5000, 0.00026, 5)

Ne_p1 <- function(Ne) {
#adjust ne for 5k
    return((Ne*5 / 5000) * 0.0066)
}
