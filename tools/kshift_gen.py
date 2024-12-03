import numpy as np
from scipy.stats import qmc

m_num = 6
num_file = 2**m_num
sampler = qmc.Sobol(d=3)
sample = sampler.random_base2(m_num)

for ifile in range(num_file):

    filename = 'kshift_coor_'+str(ifile).zfill(5)+'.out'
    with open(filename, mode='w') as f:
        out_str = 'dk_shift1 = ' + str(sample[ifile,0]) + '\n'
        f.write(out_str)
        out_str = 'dk_shift2 = ' + str(sample[ifile,1]) + '\n'
        f.write(out_str)
        out_str = 'dk_shift3 = ' + str(sample[ifile,2]) + '\n'
        f.write(out_str)
        




print(sample)


