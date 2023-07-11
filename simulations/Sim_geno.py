#!/usr/bin/env python
import argparse as ap
import sys
import numpy as np
import pandas as pd
import scipy.linalg as linalg
from pandas_plink import read_plink
from scipy import stats
import random
import warnings
import scipy
import time
import subprocess, os
import gzip
import os.path  #new
from os import path  #new

mvn = stats.multivariate_normal

def sim_geno(L, n): 
    """
    removed first argument with is L (works with regular LD matrix too), added p to indicate number of snps.
    Sample genotypes from an MVN approximation.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param n: int the number of genotypes to sample

    :return: numpy.ndarray n x p centered/scaled genotype matrix
    """
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T
    #Z = np.random.normal(size=(n, p)).T #snps x people
    #Z = Z.T #I added this. people x snps
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0) #173 x 4 in original, ine too after transpose.
    return Z

args = sys.argv
file_name = args[1] # path to plink files of simulated gene
nums_people = args[2] # number of people to simulated genotypes for 
pop = args[3] # population group for naming

print(file_name)

# Create simulated 1Kg cohort for imputation 

bim, fam, G = read_plink(file_name, verbose=False) #Note: we will want to use EUR 1KG and AFR 1KG to simulate this 
G = G.T
bim = np.array(bim)

np.random.seed(12345)    
# estimate LD for population from PLINK data
n, p = [float(x) for x in G.shape]
p_int = int(p)
mafs = np.mean(G, axis=0) / 2
G -= mafs * 2
G /= np.std(G, axis=0)
# regularize so that LD is PSD (positive semi definite)
LD = np.dot(G.T, G) / n + np.eye(p_int) * 0.1 #this is the LD matrix!!!
# compute cholesky decomp for faster sampling/simulation
L = linalg.cholesky(LD, lower=True)
Z_qtl = pd.DataFrame(sim_geno(L,int(nums_people)))  
filename = "simulated_genotypes/simulated_genotypes_" + pop + ".csv"
Z_qtl.round(4).to_csv(filename, sep="\t", index=False,mode = "a", header = False)


      
