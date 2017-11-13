#################################################
#  Implementation of BTW Sandpile simulation 	#
#  Jose Barreiros, May 2017						#
#  SYSEN6000 Cornell University					#
#################################################

import pdb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colorbar as cbar
from numpy.random import randn
from scipy import integrate
import matplotlib
from scipy import signal as sg
import subprocess                 # For issuing commands to the OS.
import numpy.matlib
import random
import csv

z_crit=4  #global variable
fig=plt.figure() #global figure

def twoD_square(L):  #create a 2D LxL matrix 
	t0=np.zeros((L,L))
	return t0

def init_pile(m):   #	randomly initialiaze matrix m from z_crit+1 to a. 
	a=9
	t1 = np.random.randint(z_crit,9+1,m.shape)   #(z_crit,3*z_crit+1,z.shape)   
	return t1

def add_boundary(p):   #add a boundary of zeros in a matrix p. Correspond to sand falling at the edge of the grid
    t2=np.lib.pad(p,1,'constant',constant_values=(0))
    return t2

def clear_boundary(q):   #clear the boundary of a  matrix. Correspond to sand falling at the edge of the grid
    size=q.shape[0]
    q[0,0:size]=np.zeros(size,dtype='int32')
    #pdb.set_trace()
    q[0:size,size-1]=np.transpose(np.matrix(np.zeros(size,dtype='int32')))
    q[0:size,0]=np.transpose(np.matrix(np.zeros(size,dtype='int32')))
    q[size-1,0:size]=np.zeros(size,dtype='int32')
    #qq=q.tolist()
    return q



def sand_rule_opt(r):   #apply rule to matrix r
	aux=r
	#print(aux)
	stencil=np.matrix([[0,1,0],[1,-4,1],[0,1,0]]) #neigbour stencil
	oness=np.greater(r,z_crit*np.ones(r.shape)).astype(numpy.int64)  #check for z[i,j]>zcrit
	
	#print(oness)
	Dz=sg.convolve(oness,stencil,'same')    #covolve stencil with the checked (z>z_crit) matrix 
	#print(Dz)
	aux+=Dz.astype(numpy.int64)

	r=clear_boundary(aux)
	#count number of sites that toppled
	#print(np.less(clear_boundary(Dz),np.zeros(Dz.shape)).astype(numpy.int64))
	#n_toppled=np.less(clear_boundary(Dz),np.zeros(Dz.shape)).astype(numpy.int64).sum()
	#
	return r

def plotting(a,j,jndex):   #plot matrix a representing j-th toppling in jndex-th avalanche

	ax=fig.add_subplot(111)
	ax.set_title('Height of the Sandpile')
	cax = ax.imshow(a, interpolation='nearest')
	cax.set_clim(vmin=0, vmax=8)
	cbar = fig.colorbar(cax, ticks=[0,3, 5, 8], orientation='vertical')
	filename = str('%01d_%03d' % (jndex,j) + '.png')
	plt.savefig(filename, dpi=100)
	print('Wrote file', filename)
	plt.clf()
	return 0


# Mencoder code copied from matplotlib website
# Mencoder stuff is copyright Josh Lifton 2004
# 'Permission is hereby granted to use and abuse this document
# so long as proper attribution is given.' Josh Lifton 2004

def movie():    #create movie with images from ploting function. 
	print("\nWriting movie\n")
	command = ('mencoder',
            'mf://*.png',
            '-mf',
            'type=png:w=800:h=600:fps=25',
            '-ovc',
            'lavc',
            '-lavcopts',
            'vcodec=mpeg4',
            '-oac',
            'copy',
            '-o',
            'output.avi')
	print("\n\nabout to execute:\n%s\n\n" % ' '.join(command))
	subprocess.check_call(command)
	print("\n\n The movie was written to 'output.avi'")

def run(s, k):  #k is index to plot
	
	i=0
	#print("\nZ[initial]:\n")
	#print(s)
	t4=np.greater(s,z_crit*np.ones_like(s))
	n_toppled=np.sum(t4)
	#print(t4)
	#print(n_toppled)
	#plotting(s,0,k)
	
	while np.array_equal(np.less_equal(s,z_crit*np.ones_like(s)),np.ones_like(s,dtype=bool))==False:  #run until z<=z_crit
		
		z2=sand_rule_opt(s)
		s=z2
		i+=1  #number of time steps to require to terminate the avalanche
		#print("\nz[%d]:\n" % (i))
		#print(s)
		t5=np.greater(s,z_crit*np.ones_like(s))
		t4+=t5
		n_toppled+=np.sum(t5)  #number of sites that toppled
		#print(t4)   #fig1. sites that toppled 
		#print(n_toppled) 
		#plotting(s,i,k)
	return i, n_toppled, 1*t4,s


def avalanche(u,index):  #Perform 1 avalanche.  index is to plot   

	zt=np.matrix(u)
	size=zt.shape[0]-2  #find the size without boundaries
	x=random.randint(0,size-1)  #calculate random i position in matrix
	y=random.randint(0,size-1)  #calculate random j position in matrix
	zt[x+1,y+1]+=15  #raise a random size over z_crit
	#print("raised Z in x:%d y:%d:\n" %(x,y))
	#print(zt)
	[time,size_avalanche,im_av,av_out]=run(zt, index)
	#plotting(im_av,-1,index)  #plot the avalanche like in Fig1
	##print('avalanche:%d\ttime to complete:%d\tsize of avalanche:%d\t\n' %(index,time,size_avalanche))
	return av_out,time,size_avalanche  #return the matrix after the avalanche

def n_avalanches(v,N):
	i=0
	zt2=v
	time_vector=[]
	size_vector=[]
	for i in range(0,N):
		
		#print("init Z  in %d:\n" %(i))
		#print(zt2)
		[t6,ti,si]=avalanche(zt2,i)
		time_vector.append(ti)
		size_vector.append(si)

		zt2=t6
		
	return time_vector,size_vector

#z=np.matrix([[3, 4, 6, 7],[8, 9, 10, 11],[12, 13, 14, 15],[16, 17, 18, 19]])
#print("z:\n")
#print(z)
#print("with stencil\n")

#print(sand_rule_opt(z)[0])
def build_freq(a_vector):  #build a frequency table:  log(count) vs log(elements) 
	t8={x:a_vector.count(x) for x in a_vector}
	elements, count = list(t8.keys()), list(t8.values())
	#print(np.log(elements))
	#print(np.log(count))
	freq_table=np.column_stack((np.log(elements),np.log(count)))
	return freq_table

def main():

	lenght=50 #size of the grid
	n=2000 #number of avalanches
	ZZ1=add_boundary(init_pile(twoD_square(lenght))) #init values in grid
	print("Z:")
	print(ZZ1)

	#tet=add_boundary([[1,2,3],[4,5,6],[0,3,9]])
	#[time,size_avalanche,im_av,av_out]=run(np.matrix(tet),0)
	#plotting(im_av,-1,1)
	#print(av_out)
	#print(avalanche(ZZ1,1))
	t_vector,s_vector=n_avalanches(ZZ1,n) 
	#print(t_vector,s_vector)
	##V_T=build_freq(t_vector)
	##V_S=build_freq(s_vector)
	#pdb.set_trace()
	V_T=np.log(t_vector)
	V_S=np.log(s_vector)

	with open('file3_b.csv', 'w') as myfile:
		wr = csv.writer(myfile, delimiter=',',quoting=csv.QUOTE_ALL)
		wr.writerow(t_vector)
		wr.writerow(s_vector)
	print('End of simulation')

	figure = plt.figure(figsize=plt.figaspect(.4),facecolor='white')
	str_t='Distribution in BTW Sandpile model\n Lifetime of Avalanche'+'\t'*25 +'Size of Avalanche\n.'
	figure.suptitle(str_t,fontsize=18)
	ax = figure.add_subplot(2,2,1)
	#ax.plot(V_T[:,0],V_T[:,1], '.', linewidth=0.1)
	ax.hist(V_T,bins='auto',range=[0, 7])
	ax.set_xlabel('log(lifetime)', fontsize=14)
	ax.set_ylabel('count', fontsize=14)
	ax = figure.add_subplot(2,2,2)
	ax.hist(V_S,bins='auto',range=[0, 10])
	#ax.plot(V_S[:,0],V_S[:,1], '.', linewidth=0.1)
	ax.set_xlabel('log(size)', fontsize=14)
	ax.set_ylabel('count', fontsize=14)


	ax = figure.add_subplot(2,2,3)
	ax.hist(t_vector,bins='auto',range=[0, 200])
	ax.set_xlabel('lifetime', fontsize=14)
	ax.set_ylabel('count', fontsize=14)

	ax = figure.add_subplot(2,2,4)
	ax.hist(s_vector,bins='auto',range=[0, 1500])
	ax.set_xlabel('size', fontsize=14)
	ax.set_ylabel('count', fontsize=14)

	plt.show()
	


	#V=build_freq(s_vector)




#run(Z)
#movie()

if __name__ == '__main__':
    main()


#Z=add_boundary(init_pile(twoD_square(m_dim)))
#i=0
#print(Z)

#print(np.less_equal(Z,z_crit*np.ones_like(Z)))
#print(np.ones_like(Z,dtype=bool))

#print(np.array_equal(np.less_equal(Z,z_crit*np.ones_like(Z)),np.ones_like(Z,dtype=bool)))




