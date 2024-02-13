#Slightly modified to suit this project:
#Basic stats, sliding window, SLIm ms output
#No divergence here

from __future__ import print_function
import libsequence
import sys
import pandas
import math
import argparse
import os
import numpy

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-masking_status', dest = 'masking_status', action='store', nargs = 1, type = str, help = 'masking status: masked/unmasked/ascertained')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, default = 100, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, default = 100, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-allele_freq_filter', dest = 'allele_freq_filter', action='store', nargs = 1, type = list, help = 'count of alleles to calculate stats for')
parser.add_argument('-binSize', dest = 'binSize', action='store', nargs = 1, default = 50, type = int, help = 'bin size in bp')#50 bp

#parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
args = parser.parse_args()
masking_status = args.masking_status[0]
chr_len =  args.regionLen[0]
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
infolder = args.input_folder[0]
outfolder = args.output_folder[0]
allele_freq_filter = args.allele_freq_filter[0]
#prefix = args.output_prefix[0]
bin_size = int(args.binSize[0])
print(type(bin_size))
#clean up allele freq list
allele_freq_filter = [int(x) for x in allele_freq_filter if x.isdigit()]

print(type(allele_freq_filter))
print(allele_freq_filter)

def calculate_mean(l_values):
    if len(l_values) > 0:
        return(str(round(numpy.mean(l_values), 3)))
    else:
        return("NA")
def refresh_dict(BIN_SIZE, TOT_LEN):
    d_LD_BINNED = {}
    d_LD_BINNED['rsq'], d_LD_BINNED['D'], d_LD_BINNED['Dprime'] = {}, {}, {}
    i = 1
    while i <= TOT_LEN:
        BIN = str(i) + "-" + str(i + BIN_SIZE - 1)
        #print(BIN)
        d_LD_BINNED['rsq'][BIN] = []
        d_LD_BINNED['D'][BIN] = []
        d_LD_BINNED['Dprime'][BIN] = []
        i = i + BIN_SIZE
    #print(d_LD_BINNED)
    return(d_LD_BINNED)

def find_my_bin(DIST, BIN_SIZE, TOT_LEN):
    if DIST == 0:
        DIST = 1 #make it so, 0 doesn't mean anything.
    #print(DIST)
    BIN = ""
    i = 1
    while i <= TOT_LEN:
        bin_start = i
        bin_end = bin_start + int(BIN_SIZE) - 1
        if i <= DIST:
            #print("yes 1")
            #print (bin_start)
            #print (bin_end)
            if int(DIST) <= bin_end:
                #print ("yes 2")
                BIN = str(bin_start) + "-" + str(bin_end)
        i = i + BIN_SIZE
    if BIN == "":
        print ("error finding bin for " + str(DIST))
    else:
        #print(BIN)
        return(BIN)

def bin_LD_stats_by_distance(l_LD_stats, BIN_SIZE, d_LD_BINNED):
    for d_LD in l_LD_stats:
        DIST = round((float(d_LD['j'])*chr_len) - (float(d_LD['i'])*chr_len))
        #print(DIST)
        d_LD_BINNED['rsq'][find_my_bin(DIST, int(BIN_SIZE), int(args.winSize[0]))].append(d_LD['rsq'])
        d_LD_BINNED['D'][find_my_bin(DIST, int(BIN_SIZE), int(args.winSize[0]))].append(d_LD['D'])
        d_LD_BINNED['Dprime'][find_my_bin(DIST, int(BIN_SIZE), int(args.winSize[0]))].append(d_LD['Dprime'])
    return()
#def get_S(f_ms):
#    for line in f_ms:
#        line1 = line.strip('\n')
#        if "segsites" in line1:
#            S = line1.split()[1]
#    f_ms.seek(0)
#    return S
#result files:
print(','.join(map(str, allele_freq_filter)))
result =  open(outfolder + "/" + "allresults_window" + "_" + str(args.winSize[0]) + "_" + ','.join(map(str, allele_freq_filter)) +"BinSize" + str(bin_size)  +  ".stats", 'w+')
result.write("filename" + '\t' + "output" + '\t' + "site" + '\t' + "posn" + '\t' + "S" + '\t' + "thetapi" + '\t' + "thetaw" + '\t' + "thetah" + '\t' + "hprime" + '\t' + "tajimasd" +  '\t' + "numSing" + '\t' + "hapdiv" + '\t' + "Concatenated" + '\t' + "rsq" + '\t' + "D" + '\t' + "Dprime" + '\n')

if os.path.exists(infolder) and os.path.isdir(infolder):
    directories = [d for d in os.listdir(infolder) if os.path.isdir(os.path.join(infolder, d))]
    print(directories)
else:
    print("The specified path is not a valid directory.")
print(directories)
for prefix in directories:
    print(prefix)

    #go through all simulation replicates and read data into pylibseq format
    #addin the option of ignoring some files if they don't exist
    #if masking_status == "masked":
    #    os.system("ls " + infolder + "/*_masked.ms > " + outfolder + "/" + prefix + ".list")
    #elif masking_status == "ascertained":
    #    os.system("ls " + infolder + "/*_singletons.ms > " + outfolder + "/" + prefix + ".list")
    #
    #f_list = open(outfolder + "/" + prefix + ".list", 'r')
    #numsim = 1
    #s_absent = 0
    #for Aline in f_list:
    #    Aline1 = Aline.strip('\n')
    #    f_name = Aline1.split(".")[0]
    #    f_name_small = Aline1.split("/").pop()
    #    print ("Reading file:" + Aline1)
    #    #try:
    #    if numsim > 0:
    for output in range(1,6):
        #print(output)
        for j in ["neutral", "selected"]:
            #print(j)
            current_output = "output" + str(output) + "_" + j
            #print(current_output)
            f_ms = open(infolder + "/" + prefix + "/" + current_output + ".ms", 'r')
            #S = get_S(f_ms)
            l_Pos = [] #list of positions of SNPs
            l_Genos = [] #list of alleles
            d_tmp = {}
            for line in f_ms:
                line1 = line.strip('\n')
                if "positions" in line1:
                    line2 = line1.split()
                    i = 0
                    for x in line2:
                        if "position" not in x:
                            l_Pos.append(float(x))
                            d_tmp[str(i)] = ""
                            i = i + 1
                elif "//" not in line and "segsites" not in line:
                    #print (d_tmp)
                    i = 0
                    while i < len(line1):
                        d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                        i = i + 1
            #print (d_tmp)
            l_data = []
            i = 0
            while i < len(l_Pos):
                l_Genos.append(d_tmp[str(i)])
                t_tmp = (l_Pos[i], d_tmp[str(i)])
                l_data.append(t_tmp)
                i = i + 1
                #print (l_Pos)
                #print (l_Genos)

            #print(l_data)
            #assign object
            sd = libsequence.SimData(l_data)
            #print(sd)
            # filter allele freqs based on provided parameter. Ex. if ==1, only analyze doubletons
            #allele_freq_filter = [1,2,3,4,5]
            if allele_freq_filter != 0:
                sd = libsequence.removeColumns(sd,lambda x : x.one in allele_freq_filter)
                #print(sd)
            #sd.assign(l_Pos[10:100],l_Genos[10:100]
            #define sliding windows:
            w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
            #chromosome length = 30kb, window size = 5 kb
            num_win = len(w)
            #set up a new dict to store LD stats from all windows:
            d_LD_binned = refresh_dict(bin_size, int(args.winSize[0]))
            #calculate summary statistic in sliding window:
            print ("calculating stats in windows")
            win_name = 1
            for i in range(len(w)):
                wi = w[i]
                #print (wi)
                pswi = libsequence.PolySIM(wi)
                result.write(prefix + '\t' + current_output + '\t' + j + '\t' + str(win_name) + '\t' + str(pswi.numpoly()) + '\t' + str(pswi.thetapi()) + '\t' + str(pswi.thetaw()) + '\t' + str(pswi.thetah()) + '\t' + str(pswi.hprime()) + '\t' + str(pswi.tajimasd()) + '\t' + str(pswi.numexternalmutations()) + '\t' + str(pswi.hapdiv()) + '\t')
                
                #read data to calculate LD based stats:
                
                #if len(wi.pos()) >= 5: #These are pairwise stats. If only 1 site exists, it'll show an error.
                    #print (i)
                LD_tmp = libsequence.ld(wi)
                bin_LD_stats_by_distance(LD_tmp, bin_size, d_LD_binned)
                win_name = win_name + 1
                #write it in a file:
                result_df = pandas.DataFrame(columns=['mean_dist', 'rsq_mean', 'D_mean', 'Dprime_mean'])
                for s_bin in d_LD_binned['rsq'].keys():
                    #print(s_bin)
                    mean_dist = (int(s_bin.split('-')[1]) + int(s_bin.split('-')[0])) / 2.0
                    rsq_mean = calculate_mean(d_LD_binned['rsq'][s_bin])
                    D_mean = calculate_mean(d_LD_binned['D'][s_bin])
                    Dprime_mean = calculate_mean(d_LD_binned['Dprime'][s_bin])

                    # Append the data to the DataFrame
                    new_row = {'mean_dist': mean_dist, 'rsq_mean': rsq_mean, 'D_mean': D_mean, 'Dprime_mean': Dprime_mean}
                    #print(new_row)
                    result_df = pandas.concat([result_df, pandas.DataFrame([new_row])], ignore_index=True)
                    
                # Calculate the mean of each column in the DataFrame
                #print(result_df)
                # Convert columns to numeric, replacing non-numeric and NaN values with NaN
                result_df = result_df.apply(pandas.to_numeric, errors='coerce')
                # Concatenate the rows of df1 into a single comma-separated string
                result_df['Concatenated'] = result_df.apply(lambda row: ','.join(map(str, row)), axis=1)
                #print(str(result_df['Concatenated']))
                # Transpose df1 and set it as columns in df2
                #df2 = df2.join(df1['Concatenated'])

                #column_means = result_df.mean()
                #print(column_means)
                #print(result_df['rsq_mean'].mean())
                # Transpose the DataFrame and write each row as a column
                transposed_df = result_df.T
                #print(transposed_df)
                # Select the row to be concatenated
                row_to_concat = transposed_df.iloc[4]  # You can choose any row by specifying its index

                # Concatenate the row with tab-separated columns
                concatenated_row = '\t'.join(map(str, row_to_concat))
                #print(concatenated_row)
                result.write(concatenated_row)
                #result.write(str(result_df['Concatenated']) + '\n')
                    #result.write(calculate_mean(d_LD_binned['rsq'][s_bin]) + '\t' + calculate_mean(d_LD_binned['D'][s_bin]) + '\t' + calculate_mean(d_LD_binned['Dprime'][s_bin]) + '\n')
                
                LDstats = pandas.DataFrame(LD_tmp)
                ##print(len(LDstats))
                if len(LDstats) > 0:
                    meanrsq = sum(LDstats['rsq'])/len(LDstats['rsq'])
                    meanD = sum(LDstats['D'])/len(LDstats['D'])
                    meanDprime = sum(LDstats['Dprime'])/len(LDstats['Dprime'])
                    result.write('\t' +str(meanrsq) + '\t' + str(meanD) + '\t' + str(meanDprime) + '\n')
                else:
                    result.write('\t' +"NA" + '\t' + "NA" + '\t' + "NA" + '\n') 
            
                #win_name = win_name + 1

                #numsim = numsim + 1

            #print ("Number of files not read:" + '\t' + str(s_absent))
            print ("Finished")