#This is to obtain the number of fixed mutations for specific types:

import sys
import argparse
import os

#parsing user given constants                                                   
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-numGenomes', dest = 'numGenomes', action='store', nargs = 1, type = int, help = 'the first x number of genomes you want to read')
parser.add_argument('-inputFolder', dest = 'inputFolder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-inputFileExt', dest = 'inputFileExt', action='store', nargs = 1, type = str, help = 'extension used for the full output files like .txt')
parser.add_argument('-outputFolder', dest = 'outputFolder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-outputPrefix', dest = 'outputPrefix', action='store', nargs = 1, type = str, help = 'full path to output file')
parser.add_argument('-mutnTypes', dest = 'mutnTypes', action='store', nargs = 1, default="m5", type = str, help = 'list of mutation types separated by only a comma')
parser.add_argument('-NAnc', dest = 'NAnc', action='store', nargs = 1, default="5000", type = int, help = 'ancestral population size')

#read input parameters
args = parser.parse_args()                                                      
in_folder = args.inputFolder[0]                                                 
input_file_ext = args.inputFileExt[0]                                           
out_folder = args.outputFolder[0]                                               
prefix = args.outputPrefix[0]                                                   
mutn_types = args.mutnTypes[0]                                                  
num_indv = args.numGenomes[0]
N_anc = args.NAnc[0]

def get_fixed_mutations(f_fixed, mutn_types, N_anc):
    s_fixed = 0
    for line in f_fixed:
        line1 = line.strip('\n')                                                
        if "#" not in line1:                                                    
            if "Mutations" not in line1:                                        
                if "m" in line1:                                                
                    line2 = line1.split()                                       
                    if len(line2) == 9:                                         
                        if line2[2] in mutn_types:
                            if int(line2[8]) >= 10*N_anc:
                                s_fixed += 1
    return(s_fixed)

#Open output file for counts
result_count = open(out_folder + "/" + prefix + "_" + str(num_indv) + "_" + mutn_types + "_count.fixed", 'w+')
result_count.write("filename" + '\t' + str(num_indv) + '\n')

#Make a list of all *ext files:
os.system("ls " + in_folder + "/*" + input_file_ext + " > " + out_folder + "/" + prefix + ".list")

f_list = open(out_folder + "/" + prefix + ".list", 'r')
for lineA in f_list:
    lineB = lineA.strip('\n')
    f_name = lineB.split("/").pop()
    print ("Reading file:" + lineB)
    f_fixed = open(in_folder + "/" + f_name, 'r')
    s_fixed = get_fixed_mutations(f_fixed, mutn_types, N_anc)
    f_fixed.close()

    #Write the full result in counts:
    result_count.write(f_name + '\t' + str(s_fixed) + '\n')

result_count.close()
f_list.close()
print("done")

    



