#2Nspq[F+(1-F)h]/(1+F)
# idk why he's using spq, but I'm assuming we can replace that with h?

find_h <- function(s,h) {
    F <- s / (2-s)
    print(F)
    #new_h <- (1+F)*h
    #new_h <- F+(1-F)*h/(1+F)
    new_h <- F+(1-F)*h
    return(new_h)
}

find_h(0, .1)
find_h(0, .25)
find_h(0, .5)
find_h(0, .75)

find_h(0.5, .1)
find_h(0.5, .25)
find_h(0.5, .5)
find_h(0.5, .75)

find_h(0.99, .1)
find_h(0.99, .25)
find_h(0.99, .5)
find_h(0.99, .75)
