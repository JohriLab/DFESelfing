
Red <- 8/12
Blue <- 5/12
RedBlue <- 1/12
Rednone <- 7/12
noneBlue <- 4/12
nonenone <- 0

#If major allele= A,B
A= Red
a  = 1-Red
B = 1-Blue
b = Blue
Ab= RedBlue
aB= nonenone
AB= Rednone
ab= noneBlue
D1= AB-A*B
print(D1)
A*(1-B)
(1-A)*(B)
#If minor allele= A,B
A= 1-Red
a  = Red
B = Blue
b = 1-Blue
Ab= nonenone
aB= RedBlue
AB= noneBlue
ab= Rednone
D1= AB-A*B
print(D1)
#If derived allele=A,B
A= 8/12
B = 5/12
AB=1/12
Ab=7/12
aB=4/12

D1= AB-A*B
print(D1)
D2 = Ab-A*B
print(D2)
D3 = aB - a*B
print(D3)
D4 = ab - a*b
print(D4)