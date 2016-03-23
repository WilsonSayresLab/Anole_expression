'''This script will perform a permutation test on the means and medians of
Ka/Ks values between two sets of values (Assumed to be Autosomal and
X-linked).

usage: permutation.py <chrx kaks file> <autosomal KaKS file>

Requires SciPy.'''

import random
from statistics import mean
from statistics import median
from sys import argv

def buildList(chrx, auto):
    # Builds separate lists fr X and autosomal Ka, Ks, and Ka,KS values
    with open(chrx, "r") as xchr:
        xka = []
        xks = []
        xkaks = []
        for line in xchr:
            # Skip header
            if line[0] == "G":
                pass
            else:
                # Append values to lists
                splt = line.split(",")
                xka.append(eval(splt[7]))
                xks.append(eval(splt[8]))
                xkaks.append(eval(splt[9]))
    with open(auto, "r") as autosomes:
        aka = []
        aks = []
        akaks = []
        for line in autosomes:
            # Skip header
            if line[0] == "G":
                pass
            else:
                # Append values to lists
                splt = line.split(",")
                aka.append(eval(splt[7]))
                aks.append(eval(splt[8]))
                akaks.append(eval(splt[9]))
    lena = len(akaks)
    meanka, medianka = permute(aka, xka, lena)
    meanks, medianks = permute(aks, xks, lena)
    meankaks, mediankaks = permute(akaks, xkaks, lena)
    print("P value for median Ka: ", medianka/10000)
    print("P value for median Ks: ", medianks/10000)
    print("P value for median Ka/Ks: ", mediankaks/10000)
    print("P value for mean Ka: ", meanka/10000)
    print("P value for mean Ks: ", meanks/10000)
    print("P value for mean Ka/Ks: ", meankaks/10000)

def permute(a, x, lena):
    meancount = 0
    mediancount = 0
    n = 0
    # Determine observed values for each pair
    obsmean = mean(x)/mean(a)
    obsmedian = median(x)/median(a)
    # Join data sets
    join = a + x
    while n < 10000:
        # Shuffle and determine permuted means and medians
        random.shuffle(join, random.random)
        pmean = (mean(join[lena + 1:]))/(mean(join[:lena + 1]))
        pmedian = (median(join[lena + 1:]))/(median(join[:lena + 1]))
        # Add 1 to counts if permuted ratio is higher than observed
        if (pmean) >= (obsmean):
            meancount += 1
        if (pmedian) >= (obsmedian):
            mediancount += 1
        n += 1
    return meancount, mediancount

def main():
    chrx = argv[1]
    auto = argv[2]
    buildList(chrx, auto)
    
if __name__ == "__main__":
    main()
