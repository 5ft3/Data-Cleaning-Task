import numpy as np
import pandas as pd

# Load and read the text files using pandas

sumStatsTraitA = pd.read_table('sumstats_trait_A.txt')
sumStatsTraitB = pd.read_table('sumstats_trait_B.txt')

# Turn the files into a pandas data frame, so we have useable columns
dfA = pd.DataFrame(sumStatsTraitA, columns = ['SNP', 'CHR','BPOS', 'A1','A2','MAF','N', 'beta_hat', 'z', 'NCHROBS', 'info'])
dfB = pd.DataFrame(sumStatsTraitB, columns = ['SNP', 'CHR','BPOS', 'A1','A2','MAF','N', 'beta_hat', 'z', 'NCHROBS', 'info'])


#function making sure all the mismatch nucleotides are changed. This function assumes everyting in column A1 is correct data and A2 
#needed to be matched to A1. The alternative decision here was to remove rows with wrong pairings, but I decided to make this assumption
#This also will remove any data from A2 that had two values such as "AA"
def matchNucleotides(dfA, dfB):   
    dfA.loc[dfA["A1"] == "A", "A2"] = "T"
    dfA.loc[dfA["A1"] == "T", "A2"] = "A"
    dfA.loc[dfA["A1"] == "G", "A2"] = "C"
    dfA.loc[dfA["A1"] == "C", "A2"] = "G"
    dfB.loc[dfB["A1"] == "A", "A2"] = "T"
    dfB.loc[dfB["A1"] == "T", "A2"] = "A"
    dfB.loc[dfB["A1"] == "G", "A2"] = "C"
    dfB.loc[dfB["A1"] == "C", "A2"] = "G"
    return dfA, dfB

#This function assumes that any MAF > 0.50000 is actually the major allel frequency and I used the inverse to give the correct MAF. 
#For example: if the data had a MAF =  0.967340, I assume this is the major frequency and the actual MAF would be = 0.032660
#if the value presented was incorrect entirely this function would not be useful
#This function also removes any negative signs
def correctMAF(dfA, dfB):
    dfA['MAF'] = dfA['MAF'].apply(lambda x : 1 - x  if x >= 0.5 else x).abs()
    dfB['MAF'] = dfB['MAF'].apply(lambda x : 1 - x  if x >= 0.5 else x).abs()

#This function founds the BPOS to drop the decimal place
def roundBPOS(dfA, dfB):
    dfA['BPOS'] = dfA['BPOS'].apply(np.round).astype('Int64')
    dfB['BPOS'] = dfB['BPOS'].apply(np.round).astype('Int64')

#This function rounds the N column since you cannot have a fraction of a person
def roundNValue(dfA, dfB):
    dfA['N'] = dfA['N'].apply(np.ceil).astype('int', errors='ignore')
    dfB['N'] = dfB['N'].apply(np.ceil).astype('int', errors='ignore')


if __name__ == '__main__':
    matchNucleotides(dfA, dfB)
    correctMAF(dfA, dfB)
    roundNValue(dfA, dfB)
    roundBPOS(dfA, dfB)
    print(dfA, dfB)

    #compares dfA(sumstats_trait_A) and dfB(sumstats_trait_B) and prints out the rows that have polymorphisms
    diff = dfA.A1.compare(dfB.A1).rename(columns={"self": "dfA", "other": "dfB"})
    print(diff)




