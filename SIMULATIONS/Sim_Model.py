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
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import time
import subprocess, os
import gzip
import os.path  
from os import path  
from sklearn.metrics import r2_score
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

    tempdir = temp + "/temp"
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
    variance = h2g/num_causal
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
    #command = 'rm {}.log'.format(out)
    #subprocess.run(command.split())
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

#compute marginal weights per snp
def marginal_weights(ge, genotypes):
    betas = []
    pvals = []

    for i in range(genotypes.shape[1]):
        x = sm.add_constant(genotypes[:, i].reshape(-1,1))
        mod = sm.OLS(ge, x)
        results = mod.fit()
        betas.append(results.params[1])
        pvals.append(results.pvalues[1])

    return betas, pvals

#take ge, genotypes and marginal betas, create top1 gene model and find r2 between predicted and actual ge
def top1(ge, genotypes, betas):
    weights = np.array(betas.copy())
    #top 1 gene model takes the highest marginally weighted snp 
    i_max = np.argmax(weights**2)
    zero = (weights != weights[i_max])
    weights[zero] = 0
    #calculate r2
    predicted_expression = np.dot(genotypes, weights)
    r2 = r2_score(ge, predicted_expression)
    #return top1 gene model weights and r2
    return weights, r2

# FUSION framework for marginal weights - alternative to marginal_weights(ge, genotypes)
def weights_marginal_FUSION(genos, pheno, beta=True):
    n_samples = len(pheno)
    
    if beta:
        eff_wgt = np.dot(genos.T, pheno) / (n_samples - 1)
    else:
        eff_wgt = np.dot(genos.T, pheno) / np.sqrt(n_samples - 1)
    
    return eff_wgt

def load_process_sumstats(file_sumstats, bim_df):
    # LOAD AND PROCESS SUMMARY STATISTICS (SUSIE POSTERIOR) FOR MAGEPRO 

    # --- LOAD IN SUM STATS
    sumstats = pd.read_csv(file_sumstats, sep = "\t")

    # --- FLIP SIGNS, SUBSET SUM STATS TO SNPS IN COMMON
    # match snps, flip signs? - unfortunately the snps present in one population isn't always present in the other
    for index, row in sumstats.iterrows():
        matched = np.where(bim_df[:,1] == row['SNP'])
        if len(matched[0]) > 0:
            if bim_df[matched, 4] != row['A1'] or bim_df[matched,5] != row['A2']:
                if bim_df[matched, 5] == row['A1'] and bim_df[matched,4] == row['A2']:
                    sumstats.loc[index, 'BETA'] = row['BETA'] * -1
                    sumstats.loc[index, 'POSTERIOR'] = row['POSTERIOR'] * -1
                else:
                    sumstats.loc[index,'BETA'] = 0
                    sumstats.loc[index, 'POSTERIOR'] = 0

    # --- get SNPs in common between bim and sumstats
    snps_bim = pd.DataFrame(bim_df[:, 1], columns=['SNP'])
    merged = sumstats.merge(snps_bim, on='SNP', how = 'inner')
    merged = merged.set_index('SNP')
    merged = merged.reindex(bim_df[:, 1])
    merged = merged.fillna(0)
    sumstats = merged.reset_index()

    wgt = np.array(sumstats['POSTERIOR'])
    pips = np.array(sumstats['PIP'])

    return wgt, pips

def top1_cv(N, genotypes, phenotype):
    #r2all also has to be reassigned - have to perform 5-fold cross validation again using top1 model
    predicted_expressions = np.zeros(N)
    kf_top1=KFold(n_splits=5)
    for train_i_top1, test_i_top1 in kf_top1.split(genotypes):
        X_train_top1, X_test_top1 = genotypes[train_i_top1], genotypes[test_i_top1]
        y_train_top1, y_test_top1 = phenotype[train_i_top1], phenotype[test_i_top1] 
        y_train_std_top1=(y_train_top1-np.mean(y_train_top1))/np.std(y_train_top1)
        y_test_std_top1=(y_test_top1-np.mean(y_test_top1))/np.std(y_test_top1)
        wgts = weights_marginal_FUSION(X_train_top1, y_train_std_top1, beta = True)
        top1_wgts, r2_top1 = top1(y_test_std_top1, X_test_top1, wgts)
        top1_wgts_train, r2_top1_train = top1(y_train_std_top1, X_train_top1, wgts)
        predicted_expressions[test_i_top1] = np.dot(X_test_top1, top1_wgts)
    lreg = LinearRegression().fit(np.array(predicted_expressions).reshape(-1, 1), phenotype)
    r2_top1_cv = lreg.score(np.array(predicted_expressions).reshape(-1, 1), phenotype) # cross validation r2
    #full model 
    wgts_full = weights_marginal_FUSION(genotypes, phenotype, beta = True)
    coef_top1_full, r2_top1_training = top1(phenotype, genotypes, wgts_full) # top1 gene models weights, total training r2 

    return r2_top1_cv, r2_top1_training, coef_top1_full


def magepro_cv(N, geno_magepro, pheno_magepro, susie_weights):
    # N = sample size 
    # geno_magepro = genotypes
    # pheno_magepro = gene expression or phenotype data 
    # susie_weights = dictionary of external datasets to use (susie posterior weights)
        # Ancestry : Array of weights
    alphas_magepro=np.logspace(-2,1.2,50,base=10)
    kf_magepro=KFold(n_splits=5)
    #store r2 from lm for each penalty
    lm_r2_allpenalty=[] 
    for penalty_magepro in alphas_magepro:
        #everyfold record predicted expression on testing set
        predicted_expressions = np.zeros(N)
        for train_index_magepro, test_index_magepro in kf_magepro.split(geno_magepro):
            y_train_magepro, y_test_magepro = pheno_magepro[train_index_magepro], pheno_magepro[test_index_magepro]
            y_train_std_magepro=(y_train_magepro-np.mean(y_train_magepro))/np.std(y_train_magepro)
            y_test_std_magepro=(y_test_magepro-np.mean(y_test_magepro))/np.std(y_test_magepro)
            X_train_magepro, X_test_magepro = geno_magepro[train_index_magepro], geno_magepro[test_index_magepro]
            #lasso for afronly weights
            lasso = lm.Lasso(fit_intercept=True)
            lasso.set_params(alpha=best_penalty,tol=1e-4,max_iter=1000)
            lasso.fit(X_train_magepro, y_train_std_magepro)
            coef_magepro = lasso.coef_
            if np.all(coef_magepro == 0):
                wgts = weights_marginal_FUSION(X_train_magepro, y_train_magepro, beta = True)
                coef_magepro, r2_top1 = top1(y_test_std_magepro, X_test_magepro, wgts)
            #prepare for ridge regression to find optimal combination of AFR gene model and EUR sumstat
            X_train_magepro2 = np.dot(X_train_magepro, coef_magepro.reshape(-1, 1))
            X_test_magepro2 = np.dot(X_test_magepro, coef_magepro.reshape(-1, 1))
            for ancestry, weights in susie_weights.items():
                X_train_magepro2 = np.hstack((X_train_magepro2, np.dot(X_train_magepro, weights.reshape(-1, 1))))
                X_test_magepro2 = np.hstack((X_test_magepro2, np.dot(X_test_magepro, weights.reshape(-1, 1))))
            ridge = Ridge(alpha=penalty_magepro) #now finding best penalty for ridge regression 
            ridge.fit(X_train_magepro2,y_train_std_magepro)
            ridge_coef = ridge.coef_
            #predict on testing
            predicted_expressions[test_index_magepro] = np.dot(X_test_magepro2, ridge_coef)
        #record r2 from lm
        lreg = LinearRegression().fit(np.array(predicted_expressions).reshape(-1, 1), pheno_magepro)
        r2_cv_penalty = lreg.score(np.array(predicted_expressions).reshape(-1, 1), pheno_magepro)
        lm_r2_allpenalty.append(r2_cv_penalty) 
    besti_magepro=np.argmax(lm_r2_allpenalty)
    bestalpha_magepro=alphas_magepro[besti_magepro]
    magepro_r2 = lm_r2_allpenalty[besti_magepro] # returning this for cv r2

    # full model for training r2 and full coef
    lasso = lm.Lasso(fit_intercept=True)
    lasso.set_params(alpha=best_penalty,tol=1e-4,max_iter=1000)
    lasso.fit(geno_magepro, pheno_magepro)
    coef_magepro = lasso.coef_
    if np.all(coef_magepro == 0):
        wgts = weights_marginal_FUSION(geno_magepro, pheno_magepro, beta = True)
        coef_magepro, r2_top1 = top1(pheno_magepro, geno_magepro, wgts)
    X_magepro = np.dot(geno_magepro, coef_magepro.reshape(-1, 1))
    for ancestry, weights in susie_weights.items():
        X_magepro = np.hstack((X_magepro, np.dot(geno_magepro, weights.reshape(-1, 1))))
    ridge = Ridge(alpha=bestalpha_magepro)
    ridge.fit(X_magepro, pheno_magepro)
    wgts_sep = coef_magepro.reshape(-1, 1)
    for ancestry, weights in susie_weights.items():
        wgts_sep = np.hstack((wgts_sep, weights.reshape(-1, 1)))
    magepro_coef = np.dot(wgts_sep, ridge.coef_) #magepro coef

    return magepro_r2, magepro_coef

args = sys.argv
plink_file = args[1] # path to plink files of simulated gene
samplesizes = int(args[2]) # number of people
pop = args[3] #population
genotype_file = args[4] #path to simulated genotypes
sumstats_file = args[5] #path to sumstats created in other poputation
sumstats_files = sumstats_file.split(',')
population_sumstat = args[6]
populations_sumstat = population_sumstat.split(',') # comma separated list of ancestries of sumstats
set_num_causal = int(args[7]) #number of causal snps
CAUSAL = int(args[8]) # index of causal snp
sim = int(args[9]) #iteration of simulation
set_h2 = float(args[10]) #predetermiend h2g
temp_dir = args[11] #temporary directory for gcta
out_results = args[12]

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
#make sure to pick a 1 mb region with enough snps (maybe > 30?) up to you! 
betas = sim_effect_sizes(set_h2, set_num_causal, bim.shape[0], CAUSAL)
beta_causal = betas[CAUSAL]
b_qtls = np.array(betas)
b_qtls = np.reshape(b_qtls.T, (bim.shape[0], 1))
#print(b_qtls.shape)

# --- SIMULATE LASSO GENE MODEL
best_penalty, coef, h2g, hsq_p, r2all, gexpr = sim_eqtl(z_eqtl, samplesizes, b_qtls, float(set_h2), temp_dir) #this runs gcta and lasso 
#gexpr already standardized

# --- if the best lasso model with the best penalty gives coef of all 0, we have to use top1 to compute r2 afr
if np.all(coef == 0):
    print("using top1 backup")
    r2all, r2_top1, coef = top1_cv(samplesizes, z_eqtl, gexpr)

print(("heritability " + pop + ": "))
print(h2g)

# --- PROCESS SUMMARY STATISTICS
num_causal_susie = 0
sumstats_weights = {}
for i in range(0,len(sumstats_files)):
    varname = populations_sumstat[i]+'weights'
    weights, pips = load_process_sumstats(sumstats_files[i], bim)
    sumstats_weights[varname] = weights
    if pips[CAUSAL] >= 0.95:
        num_causal_susie = num_causal_susie + 1


# --- RUN CV MAGEPRO 
magepro_r2, magepro_coef = magepro_cv(samplesizes, z_eqtl, gexpr, sumstats_weights)

print("magepro: ")
print(magepro_r2)
print("afronly: ")
print(r2all)

# --- POPULATE OUTPUTS
#coef = afr lasso model gene model 
lasso_causal_nonzero = 0 # 0 for causal has 0 weight in lasso gene model. 1 for nonzero
if coef[CAUSAL] != 0:
    lasso_causal_nonzero = 1
#magepro_coef = magepro gene model 
magepro_causal_nonzero = 0 # 0 for causal has 0 weight in magepro gene model. 1 for nonzero
if magepro_coef[CAUSAL] != 0:
    magepro_causal_nonzero = 1
#beta of causal snp from afr lasso
afr_B_causal = (coef[CAUSAL])
#beta of causal snp from magepro
magepro_B_causal = (magepro_coef[CAUSAL])
#actual beta
true_B_causal = beta_causal

#h2g = afr h2
#magepro_r2 = magepro cv r2 
#r2all = afronly magepro cv r2

#filename = "results/magepro_results_" + str(samplesizes) + "_h" + str(set_h2) + ".csv"
filename = out_results + "/magepro_results_" + str(samplesizes) + "_h" + str(set_h2) + ".csv"

output = pd.DataFrame({'sim': sim, 'afr_h2': h2g, 'lasso_causal': lasso_causal_nonzero, 'magepro_causal': magepro_causal_nonzero, 'afr_beta_causal': afr_B_causal, 'magepro_beta_causal': magepro_B_causal, 'true_B_causal': true_B_causal, 'afr_r2': r2all, 'magepro_r2': magepro_r2, 'causal_susie': num_causal_susie}, index=[0])
if sim == 1:
    output.to_csv(filename, sep="\t", index=False, header = True)
else:
    output.to_csv(filename, sep="\t", index=False, header = False, mode='a')

