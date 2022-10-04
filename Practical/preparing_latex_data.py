# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 09:39:08 2022

@author: admin
"""
import scipy.io
import numpy as np
import matplotlib.pyplot as plt

data_file = scipy.io.loadmat('D://SIFA/blind sefa/Data1/HD7.mat')
Key_rank_ineff=data_file['Key_rank_ineff']
# Key_rank_joint=data_file['Key_rank_joint']
Key_rank_joint_HD=data_file['Key_rank_joint_HD']
no_traces=data_file['no_traces']
s1=(np.mean(Key_rank_ineff,1))
# s2=np.mean(Key_rank_joint,1)
s3=np.mean(Key_rank_joint_HD,1)


sampling_rate=10
no_temp=int(no_traces[0][0]/sampling_rate)
Sifa=np.zeros(no_temp)
Sefa=np.zeros(no_temp)
n_sample=np.zeros(no_temp)
n=0
for i in range (no_temp):
    Sifa[i]=s1[n]
    Sefa[i]=s3[n]
    n_sample[i]=n
    n=n+sampling_rate
combined=np.zeros((no_temp,3))
combined[:,0]=n_sample
combined[:,1]=Sifa
combined[:,2]=Sefa

np.savetxt("combined7.csv", combined, delimiter=",")


# plt.plot(Sifa)
# # plt.plot(s2)
# plt.plot(Sefa)

# print(np.mean(Key_rank_joint,1))