import numpy as np
from scipy import stats
from scipy.stats import pearsonr
mvn = stats.multivariate_normal

def random_correlated_effect(effect1, heritability, correlation):
    # effect1 = first effect for which you want another effect that is correlated to it
    # heritability = preset heritability value (variance of the distribution of effect size)
    # correlation = desired correlation of the two effect sizes 
    effect_temp = np.random.normal(0, np.sqrt(heritability))
    effect2 = (correlation*effect1) + np.sqrt(1-(correlation**2))*effect_temp # see supplements for more information
    return (effect2)

heritability = 0.05
correlation = 0.95
iteratioons = 10000
effects1 = []
effects2 = []
for i in range(iteratioons):
    effect1 = np.random.normal(0, np.sqrt(heritability))
    effects1.append(effect1)
    effects2.append(random_correlated_effect(effect1, heritability, correlation))

correlation_coefficient, _ = pearsonr(effects1, effects2)
variance = np.var(effects2)
print("Correlation between two effects: " + str(correlation_coefficient))
print("Variance of computed effect: " + str(variance))