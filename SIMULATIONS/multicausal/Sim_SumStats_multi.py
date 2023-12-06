#!/usr/bin/env python
import argparse as ap
import sys
import numpy as np
import pandas as pd
import scipy.linalg as linalg
from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm
import random
import warnings
import scipy
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.exceptions import ConvergenceWarning
import time
import subprocess, os
import gzip
import os.path  
from os import path  
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression

mvn = stats.multivariate_normal

# --- SIMULATED GENE EXPRESSION
def sim_trait(g, h2g): 
    """
    Simulate a complex trait as a function of latent genetic value and env noise.

    :param g: numpy.ndarray of latent genetic values
    :param h2g: float the heritability of the trait in the population

    :return: (numpy.ndarray, float) simulated phenotype, sd of Y
    """
    n = len(g)
    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e #appears to add noise to the genetic values, adding the same about of noise if the seed is the same 
    else:
        e = np.random.normal(0, 1, n)
        y = e
    # standardize
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()
    return y, y_std  

# --- FIT GENE MODEL
def fit_lasso_katie(Z, y):
    """
    Infer eqtl coefficients using LASSO regression for simeQTL function above. Uses the PLINK-style coordinate descent algorithm
    that is bootstrapped by the current h2g estimate.
    :param Z:  numpy.ndarray n x p genotype matrix
    :param y: numpy.ndarray gene expression for n individuals
    :return: (numpy.ndarray, float, float) tuple of the LASSO coefficients, the r-squared score, and log-likelihood
    """
    n, p = Z.shape
    alphas=np.logspace(-2,1.2,50,base=10)
    kf=KFold(n_splits=5)

    #lm() r2 extraction
    lm_r2_allpenalty=[]

    for penalty in alphas:

        #every fold, predict expression
        predicted_expressions = np.zeros(n)

        CV_r2=[]
        for train_index, test_index in kf.split(Z):
            X_train, X_test = Z[train_index], Z[test_index]
            y_train, y_test = y[train_index], y[test_index]
            y_train_std=(y_train-np.mean(y_train))/np.std(y_train)
            y_test_std=(y_test-np.mean(y_test))/np.std(y_test)
            lasso = lm.Lasso(fit_intercept=True)
            lasso.set_params(alpha=penalty,tol=1e-4,max_iter=1000)
            lasso.fit(X_train,y_train_std)

            #predict on testing
            predicted_expressions[test_index] = np.dot(X_test, lasso.coef_)

        #evaluate r2 using lm
        lreg = LinearRegression().fit(np.array(predicted_expressions).reshape(-1, 1), y)
        r2_cv_penalty = lreg.score(np.array(predicted_expressions).reshape(-1, 1), y)
        lm_r2_allpenalty.append(r2_cv_penalty)

    #print(CV_r2_allpenalty)
    #besti=np.argmax(CV_r2_allpenalty)
    besti=np.argmax(lm_r2_allpenalty)
    bestalpha=alphas[besti]
    lasso = lm.Lasso(fit_intercept=True)
    lasso.set_params(alpha=bestalpha,tol=1e-4,max_iter=1000)
    lasso.fit(Z, y)
    coef = lasso.coef_
    return bestalpha, coef, lm_r2_allpenalty[besti]

# --- SIMULATE EQTL
def sim_eqtl(Z_qtl, nqtl, b_qtls, eqtl_h2, temp):
    """
    Simulate an eQTL study using `nqtl` individuals.
    :param Z_qtl: simulated genotypes
    :param nqtl: int the number of eQTL-panel genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene: SNPs x genes
    :param eqtl_h2: list of expression variance explained by linear model of SNPs for each gene (not used) rather min(1,b@LD_qtl[np.ix_(nearSNPsi,nearSNPsi)]@b)
    :param temp: name of temporary directory
    :return:  (numpy.ndarray x 6) Arrays of best lasso penalty, eqtl effect sizes, heritability, cross validation accuracies, simulated ge
    """
    allqtlshape = b_qtls.shape
    n, p = [float(x) for x in Z_qtl.shape]

    tempdir = "/expanse/lustre/projects/ddp412/kakamatsu/eQTLsummary/multipopGE/simulations/" + temp + "/temp"
    b = b_qtls[:,0]
    nonzeroindex = np.where(b != 0)[0]
    gexpr = sim_trait(np.dot(Z_qtl,b), eqtl_h2)[0]
    penalty, coef, r2= fit_lasso_katie(Z_qtl, gexpr)
    if r2>0:
        h2g,hsq_se,hsq_p = get_hsq(Z_qtl,gexpr,tempdir,"/expanse/lustre/projects/ddp412/kakamatsu/fusion_twas-master/gcta_nr_robust")
    else:
        h2g=0
        hsq_se=10
        hsq_p=1
    return (penalty, coef,h2g,hsq_p,r2,gexpr)

# -- SIMULATE CAUSAL EFFECT SIZES
def sim_effect_sizes(h2g, num_causal, num_snps, causal_index):
    mean = 0
    variance = h2g/np.sqrt(num_causal)
    effect_sizes = np.zeros(num_snps)
    if isinstance(causal_index, list):
        for index in causal_index:
            effect_sizes[index] = np.random.normal(mean, np.sqrt(variance))
    else:
        effect_sizes[causal_index] = np.random.normal(mean, np.sqrt(variance))
    return effect_sizes

def get_hsq(geno, gene_expr, out, gcta):
    nindv, nsnp = geno.shape
    grm = np.dot(geno, geno.T) / nsnp
    write_pheno(gene_expr, out)
    write_grm(grm, nsnp, out)
    run_gcta(gcta, out)
    clean_up(out)
    hsq_f = '{}.hsq'.format(out)
    if os.path.exists(hsq_f):
        hsq_out = pd.read_table(hsq_f)
        hsq = hsq_out['Variance'][0]
        hsq_se = hsq_out['SE'][0]
        hsq_p = hsq_out['Variance'][8]
        command = 'rm {}'.format(hsq_f)
        subprocess.run(command.split())
        return hsq, hsq_se, hsq_p
    else:
        return 0,10000,1
    return None

def run_gcta(gcta, out):
    command = '{} --grm-gz {} --pheno {} --reml --reml-no-constrain --thread-num 4 --out {}'\
        .format(gcta, out, out+'.phen', out)
    subprocess.run(command.split())
    command = 'rm {}.phen'.format(out)
    subprocess.run(command.split())

def write_pheno(phen, out):
    out_f = open(out+'.phen', 'w')
    nindv = phen.shape[0]
    for i in range(nindv):
        out_f.write('{}\t{}\t{}\n'.format(i, i, phen[i]))
    out_f.close()

def clean_up(out):
    command = 'rm {}.grm.gz'.format(out)
    subprocess.run(command.split())
    command = 'rm {}.grm.id'.format(out)
    subprocess.run(command.split())

def write_grm(grm, nsnp, out):
    out_id_f = open(out+'.grm.id', 'w')
    out_f = gzip.open(out+'.grm.gz', 'w')
    nindv = grm.shape[0]
    for i in range(nindv):
        for j in range(0, i+1):
            out_f.write('{}\t{}\t{}\t{}\n'.format(i+1, j+1, nsnp,
                grm[i,j]).encode('utf-8'))
        out_id_f.write('{}\t{}\n'.format(i, i))
    out_f.close()
    out_id_f.close()


args = sys.argv
plink_file = args[1] # path to plink files of simulated gene
samplesizes = int(args[2]) # number of people to simulated genotypes for 
pop = args[3] #population
genotype_file = args[4] #path to simulated genotypes
set_num_causal = int(args[5]) #how many causal variants?
CAUSAL = [ int(index) for index in args[6].split(',') ]#indices of causal variants
sim = int(args[7]) #iteration of simulation
set_h2 = float(args[8]) #predetermined h2g
samplesize_target = int(args[9]) #sample size
out_sumstat = args[10] #output dir for sumstats
temp_dir = args[11] #dir for temporary files for gcta

# --- READ IN SIMULATED GENE PLINK FILES
bim, fam, G_all = read_plink(plink_file, verbose=False)  
G_all = G_all.T
bim = np.array(bim)

# --- READ IN SIMULATED GENOTYPE
z_eqtl = np.array(pd.read_csv(genotype_file, sep = "\t", header = None))

# --- SIMULATE CAUSAL EFFECT SIZES 
#causal eqtl effect sizes; matrix dim = snps x 1 gene
#if our sim only has 1 gene, then we just need a simple vector of mostly zeros and 1 nonzero effect size if 1 causal variant. 
#draw from mean 0 var = 0.05 if 1 causal variant, draw from mean 0 var = h2g/sqrt(#causal sps) if more than 1 causal variant. 
betas = sim_effect_sizes(set_h2, set_num_causal, bim.shape[0], CAUSAL)
b_qtls = np.array(betas) 
b_qtls = np.reshape(b_qtls.T, (bim.shape[0], 1))

best_penalty, coef, h2g, hsq_p, r2all, gexpr = sim_eqtl(z_eqtl, samplesizes, b_qtls, float(set_h2), temp_dir) #this runs gcta and lasso

#check heritability
print(("heritability " + pop + ": "))
print(h2g)


# --- RUN LINEAR REGRESSION PER SNP TO CREATE SUMSTATS
ge_regression = gexpr #already standardized

betas = []
pvals = []

for i in range(z_eqtl.shape[1]):
    x = sm.add_constant(z_eqtl[:, i].reshape(-1,1))
    mod = sm.OLS(ge_regression, x)
    results = mod.fit()
    betas.append(results.params[1])
    pvals.append(results.pvalues[1])

sumstats = pd.DataFrame({'snp':bim[:,1], 'beta':betas, 'p':pvals, 'effect_A': bim[:,4], 'alt_A': bim[:,5]})

filename = out_sumstat + "/sumstats_" + pop + ".csv"
sumstats.to_csv(filename, sep="\t", index=False, header = True)

lasso_causal_nonzero = len([index for index in CAUSAL if coef[index] != 0]) #number of nonzero causal variants 

#write results
filename = "results_multicausal/eur_h2_causal_" + str(samplesize_target) + "_h" + str(set_h2) + ".csv"
eur_h2_causal = pd.DataFrame({'sim': sim, 'eur_h2': h2g, 'eur_lasso_casual': lasso_causal_nonzero, 'eur_r2': r2all}, index=[0])
if sim == 1:
    eur_h2_causal.to_csv(filename, sep="\t", index=False, header = True)
else:
    eur_h2_causal.to_csv(filename, sep="\t", index=False, header = False, mode='a')

