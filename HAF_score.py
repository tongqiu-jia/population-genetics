import numpy as np

with open('sample_data/haf.txt') as f:
    m = []
    for line in f:
        line = line.strip()
        row = [int(c) for c in line]
        m.append(row)
f.close()
mat = np.asarray(m)


'''
This is a simple implementation of the HAF score discribed in 
	Ronen R, Tesler G, Akbari A, Zakov S, Rosenberg NA, Bafna V. 
	Predicting Carriers of Ongoing Selective Sweeps without Knowledge of the Favored Allele. 
	Coop G, ed. PLoS Genetics. 2015;11(9):e1005527. doi:10.1371/journal.pgen.1005527.
'''
def HAF_score(snp):
    w_all = np.sum(snp,axis=0)
    
    #find haplotype
    HAF_i = []
    HAF = []
    for i in range(snp.shape[0]):
        for j in range(snp.shape[1]):
            if snp[i][snp.shape[1]-1-j] == 1:
                if (snp.shape[1]-1-j) not in HAF_i:
                    HAF_i.append((i,snp.shape[1]-1-j))
                    break
    #HAF score for each of the haplotype
    for leaf in HAF_i:
        h = snp[leaf[0]]
        c = np.multiply(w_all,h)
        HAF.append(np.sum(c))
    return HAF
