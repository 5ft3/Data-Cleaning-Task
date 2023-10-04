import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
import statsmodels.api as sm
from statsmodels.distributions.empirical_distribution import ECDF
import math


sumStatsTraitA = pd.read_table('sumstats_trait_A.txt')
sumStatsTraitB = pd.read_table('sumstats_trait_B.txt')

# Turn the files into a pandas data frame, so we have useable columns
dfA = pd.DataFrame(sumStatsTraitA, columns = ['SNP', 'CHR','BPOS', 'A1','A2','MAF','N', 'beta_hat', 'z', 'NCHROBS', 'info'])
dfB = pd.DataFrame(sumStatsTraitB, columns = ['SNP', 'CHR','BPOS', 'A1','A2','MAF','N', 'beta_hat', 'z', 'NCHROBS', 'info'])

#Calculating and printing p-values for sumstats_A and sumstats_B

p_value_dfa = scipy.stats.norm.sf(abs(dfA["z"]))*2 

# scaling down pvalues to -log10 (- log base 10 (x) = log base 10 (1/x))
log_x_scale_dfa = [math.log(1/x_i) for x_i in p_value_dfa]

# get empirical data via library function
empirical_CDF_dfa = ECDF(log_x_scale_dfa)

# plot empirical CDF data x and y in a scatter
plt.scatter(empirical_CDF_dfa.x, empirical_CDF_dfa.y)
plt.legend(["sumStatsTraitA"])

# repeat for study 2 (dfB)
p_value_dfb = scipy.stats.norm.sf(abs(dfB["z"]))*2 
log_x_scale_dfb = [ math.log(1/x_i) for x_i in p_value_dfb]
empirical_CDF_dfb = ECDF(log_x_scale_dfb)

# plot again, tried to put a legend for dfB but ran out of time
plt.scatter(empirical_CDF_dfb.x, empirical_CDF_dfb.y)
plt.legend(["sumStatsTraitB"])

# scale down the graph view
plt.xscale('log')
plt.yscale('log')

# rename labels
plt.xlabel('p-values under the null hypothesis')
plt.ylabel('empirical CDF of the p-values from the GWAS summary statistics')

# build 45 degree solid line
xpoints = ypoints = plt.xlim()
plt.plot(xpoints, ypoints, linestyle='-', color='k', lw=1)

# show graph
plt.show()


