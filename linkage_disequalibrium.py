import numpy as np
from scipy import stats
import pylab as plt

#read SNP matrix
f = open("sample.txt",'rb')
m1 = []
for line in f:
    line = line.strip()
    row = [int(c) for c in line]
    m1.append(row)
f.close()
snp = np.asarray(m1)

'''
Linkage (dis)equalibrium:
Given two bi_allelic sites with allele 0 and allele 1, 
Define P00 = Pr(allele 0 in locus 1 and locus 2)
	   P0_ = Pr(allele 0 in locus 1)
	   P_0 = Pr(allele 0 in locus 2)
We say there is linkage equalibrium if P00 = P0_ * P_0
The D-measurement of LD: 
	D = P00 - P0_ * P_0
The D' (D_prime) measurement of LD is obtained by the divding D 
	by the largest possible value:
	D' = D/Dmax
'''
def D_prime(array):
    #initialization
    num_row = array.shape[0]
    num_col = array.shape[1]
    c00 = np.zeros([num_col,num_col])
    c0_ = np.zeros([num_col,num_col])
    c_0 = np.zeros([num_col,num_col])
    c1_ = np.zeros([num_col,num_col])
    c_1 = np.zeros([num_col,num_col])
    D_prime = np.zeros([num_col,num_col])

    #find count
    for i in range(num_col):
        for j in range(num_col):
            curr_pair = array[:,[i,j]]
            for pair in curr_pair:
                if(pair[0] == 0 and pair[1] == 0):
                    c00[i][j] += 1
                if(pair[0] == 0):
                    c0_[i][j] += 1
                if(pair[0] == 1):
                    c1_[i][j] += 1
                if(pair[1] == 0):
                    c_0[i][j] += 1
                if(pair[1] == 1):
                    c_1[i][j] += 1
    #calculate probabilities: counts divided num of individuals
    p00 = c00 / num_row
    p0_ = c0_ / num_row
    p_0 = c_0 / num_row
    p1_ = c1_ / num_row
    p_1 = c_1 / num_row

    #LD
    D = np.subtract(p00,np.multiply(p0_,p_0))

    #Normalization: D' measurement
    for i in range(D.shape[0]):
        for j in range(D.shape[1]):
            if D[i][j] > 0:
                D_prime[i][j] = D[i][j]*1.0 / min(p0_[i][j]*p_1[i][j],p1_[i][j]*p_0[i][j]) 
            if D[i][j] < 0:
                D_prime[i][j] = D[i][j]*1.0 / max(-(p0_[i][j]*p_0[i][j]),-(p1_[i][j]*p_1[i][j]))
    
    # Calculation of p-value
    r = np.zeros([num_col,num_col])
    r = D / np.power(np.multiply(np.multiply(p1_,p0_),np.multiply(p_1,p_0)),1/2)
    chi_square = np.multiply(np.power(r,2),num_col)
    p_val = stats.chi2.sf(chi_square,1)

    D_alt = -np.log(p_val)
    return D_prime, D_alt


if __name__ == '__main__':
	d_prime,d_alt = D_prime(snp)

	#Show heapmap of DL using D-prime measure 
	plt.imshow(d_prime)
	plt.colorbar() 
	plt.show()