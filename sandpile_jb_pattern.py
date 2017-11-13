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

z_crit=4  #global variable
fig=plt.figure() #global figure

def twoD_square(L):  #create a 2D LxL matrix 
	t0=np.zeros((L,L))
	return t0

def init_pile(m):   #	randomly initialiaze matrix m from z_crit+1 to a*z_crit, where a=3. 
	t1 = np.random.randint(z_crit,9+1,m.shape)   #(z_crit,3*z_crit+1,z.shape)   
	return t1

def add_boundary(p):   #write a boundary of zeros in a matrix p. Correspond to sand falling at the edge of the grid
    t2=np.lib.pad(p,1,'constant',constant_values=(0))
    return t2

def clear_boundary(q):   #write a boundary of zeros in a matrix. Correspond to sand falling at the edge of the grid
    size=q.shape[0]
    q[0,0:size]=np.zeros(size,dtype='int32')
    #pdb.set_trace()
    q[0:size,size-1]=np.transpose(np.matrix(np.zeros(size,dtype='int32')))
    q[0:size,0]=np.transpose(np.matrix(np.zeros(size,dtype='int32')))
    q[size-1,0:size]=np.zeros(size,dtype='int32')
    #qq=q.tolist()
    return q

##def sand_rule(r):   #apply rule to matrix r
##	table_size=z.shape[0]
##	z2=z.copy()
##	z=np.lib.pad(z2,1,'constant',constant_values=(0))  #padded with zeros in border
##	z2=z.copy()
##	print(z_crit)
	
##	for x in range(1,table_size+1): 
##		for y in range(1,table_size+1):
			#pdb.set_trace()
			#if z2[x,y] > z_crit:  #check for colapse
##			if 2>1:
##				z[x,y]-=4   #colapse
##				z[x+1,y]+=1
##				z[x-1,y]+=1
##				z[x,y+1]+=1
##				z[x,y-1]+=1
##	return z

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

def movie():
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

def run(s, k):  #k to plot
	
	i=0
	print("\nZ[initial]:\n")
	print(s)
	t4=np.greater(s,z_crit*np.ones_like(s))
	n_toppled=np.sum(t4)
	print(t4)
	print(n_toppled)
	plotting(s,0,k)
	#print(np.array_equal(np.less_equal(z,z_crit*np.ones_like(z)),np.ones_like(z,dtype=bool)))
	while np.array_equal(np.less_equal(s,z_crit*np.ones_like(s)),np.ones_like(s,dtype=bool))==False:  #run until z<=z_crit
		
		z2=sand_rule_opt(s)
		s=z2
		i+=1  #number of time steps to require to terminate the avalanche
		print("\nz[%d]:\n" % (i))
		print(s)
		t5=np.greater(s,z_crit*np.ones_like(s))
		t4+=t5
		n_toppled+=np.sum(t5)  #number of sites that toppled
		print(t4)   #fig1. sites that toppled 
		print(n_toppled) 
		plotting(s,i,k)
	return i, n_toppled, 1*t4,s


def avalanche(u,index):  #Perform 1 avalanche.  index is to plot   
	#print("\nZ[initial]:\n")
	#print(z1)
	zt=np.matrix(u)
	size=zt.shape[0]-2  #find the size without boundaries
	x=random.randint(0,size-1)  #calculate random i position in matrix
	y=random.randint(0,size-1)  #calculate random j position in matrix
	zt[x+1,y+1]=15  #raise a random size over z_crit
	print("raised Z in x:%d y:%d:\n" %(x,y))
	print(zt)
	[time,size_avalanche,im_av,av_out]=run(zt, index)
	plotting(im_av,-1,index)  #plot the avalanche like in Fig1
	print('avalanche:%d\ttime to complete:%d\tsize of avalanche:%d\t\n' %(index,time,size_avalanche))
	return av_out,time,size_avalanche  #return the matrix after the avalanche

def n_avalanches(v,N):
	i=0
	zt2=v
	time_vector=[]
	size_vector=[]
	for i in range(0,N):
		
		print("init Z  in %d:\n" %(i))
		print(zt2)
		[t6,ti,si]=avalanche(zt2,i)
		time_vector.append(ti)
		size_vector.append(si)
		print(t6)
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
	print(np.log(elements))
	print(np.log(count))
	freq_table=np.column_stack((np.log(elements),np.log(count)))
	return freq_table

def main():

	lenght=5 #size of the grid
	n=10 #number of avalancehs
	ZZ1=add_boundary(init_pile(twoD_square(lenght))) #init values in grid
	print("Z:")
	print(ZZ1)

	#tet=add_boundary([[1,2,3],[4,5,6],[0,3,9]])
	#[time,size_avalanche,im_av,av_out]=run(np.matrix(tet),0)
	#plotting(im_av,-1,1)
	#print(av_out)
	#print(avalanche(ZZ1,1))
	t_vector,s_vector=n_avalanches(ZZ1,n) 
	print(t_vector,s_vector)
	V_T=build_freq(t_vector)
	V_S=build_freq(s_vector)
	#pdb.set_trace()



	figure = plt.figure(figsize=plt.figaspect(.4),facecolor='white')
	figure.suptitle('Distribution in BTW Sandpile model\n Lifetime of Avalanche\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tSize of Avalanche\n',fontsize=18)
	ax = figure.add_subplot(1,2,1)
	ax.plot(V_T[:,0],V_T[:,1], '.', linewidth=0.1)
	ax.set_xlabel('lifetime', fontsize=14)
	ax.set_ylabel('count', fontsize=14)
	ax = figure.add_subplot(1,2,2)
	ax.plot(V_S[:,0],V_S[:,1], '.', linewidth=0.1)
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




