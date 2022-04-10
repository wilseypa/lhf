import numpy as np
import persim
import matplotlib.pyplot as plt
import pandas as pd
import os
Original = pd.read_csv ('outputPolytopal.csv')
Original.columns =['Dimension', 'Birth', 'Death']
x = os.listdir('Delaunay/')
for file1 in x:
	filej = "Delaunay/"+file1
	testing = pd.read_csv (filej)
	testing.columns =['Dimension', 'Birth', 'Death']
	for dim in [0,1,2]:
				dimVRdgm = Original[Original["Dimension"]==dim]
				dimbetadgm = testing[testing["Dimension"]==dim]
				dimVRdgm = dimVRdgm[dimVRdgm['Death']<2]
				dimbetadgm = dimbetadgm[dimbetadgm['Death']<2]
				dgm1 = np.array(np.transpose([dimVRdgm['Birth'],dimVRdgm['Death']]))
				dgm2 = np.array(np.transpose([dimbetadgm['Birth'],dimbetadgm['Death']]))
				print("dimension",dim,file1,persim.sliced_wasserstein(dgm1,dgm2),dimVRdgm.size/3,dimbetadgm.size/3)
