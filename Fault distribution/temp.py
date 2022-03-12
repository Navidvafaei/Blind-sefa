# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from aescode import dec_to_HW
from aescode import subword

import numpy as np
import matplotlib.pyplot as plt
fd = np.random.uniform(0.0,0.4,size=256)
fd = fd/np.sum(fd)
# plt.plot(fd)
s=np.zeros((256,256))
temp=np.zeros((9,256))
Pr=np.zeros((9,256))
for key in range(256):
    for i in range (256):
        s[i,key]=dec_to_HW ((subword(i^key)))
        temp[int(s[i,key]),key]=temp[int(s[i,key]),key]+1
        Pr[int(s[i,key]),key]=fd[i]+Pr[int(s[i,key]),key]
plt.plot(Pr[:,1:250])















