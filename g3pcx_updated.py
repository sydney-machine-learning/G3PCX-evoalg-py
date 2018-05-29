# !/usr/bin/python

 
# Rohitash Chandra, Centre for Translational Data Science
# University of Sydey, Sydney NSW, Australia.  2017 c.rohitash@gmail.conm
# https://www.researchgate.net/profile/Rohitash_Chandra


import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(suppress=True) 
import random
import time
import math
import os
import shutil




 
class evolution:
	def __init__(self, pop_size, dimen, max_evals,  max_limits, min_limits):
		
		self.EPSILON = 1e-20  # convergence
		self.sigma_eta = 0.01
		self.sigma_zeta = 0.01

		self.children = 2
		self.num_parents = 3 
		self.family = 2
		self.sp_size = self.children + self.family 



		self.population =   np.random.randn( pop_size  , dimen)  * 5  #[SpeciesPopulation(dimen) for count in xrange(pop_size)]
		self.sub_pop =  np.random.randn( self.sp_size , dimen )  *5  #[SpeciesPopulation(dimen) for count in xrange(NPSize)]

		self.fitness = np.random.randn( pop_size)
		self.sp_fit  = np.random.randn(self.sp_size)
		
		self.best_index = 0
		self.best_fit = 0
		self.worst_index = 0
		self.worst_fit = 0

		self.rand_parents =  self.num_parents 

		self.temp_index =  np.arange(0, pop_size)

		self.rank =  np.arange(0, pop_size)

		self.list = np.arange(0, self.sp_size)
		self.parents = np.arange(0, pop_size)
		self.pop_size = pop_size
		self.dimen = dimen 

		self.num_evals = 0

		self.max_evals = max_evals

		self.problem = 1



 

	def fit_func(self, x):    #   function  (can be any other function, model or even a neural network)
		fit = 0

		if self.problem == 1: # rosenbrock
			for j in range(x.size -1): 
				fit  = fit +  (100.0*(x[j]*x[j] - x[j+1])*(x[j]*x[j] - x[j+1]) + (x[j]-1.0)*(x[j]-1.0))

		elif self.problem ==2:  # ellipsoidal - sphere function
			for j in range(x.size):
				fit= fit + ((j+1)*(x[j]*x[j]))
 
				  

		return fit # note we will maximize fitness, hence minimize error 

	def rand_normal(self, mean, stddev):

		n2 = 0.0 
		n2_cached = False 
		
		if (not n2_cached):
			#choose a point x,y in the unit circle uniformly at random

			x = np.random.uniform(-1,1,1) 
			y = np.random.uniform(-1,1,1)  
			r = x*x + y*y 

			while (r == 0 or r > 1):
				x = np.random.uniform(-1,1,1) 
				y = np.random.uniform(-1,1,1)  
				r = x*x + y*y 

			# Apply Box-Muller transform on x, y
			d = np.sqrt(-2.0*np.log(r)/r) 
			n1 = x*d
			n2 = y*d
			# scale and translate to get desired mean and standard deviation

			result = n1*stddev + mean  
			n2_cached = True 
			return result


		else:
			n2_cached = False
			return n2*stddev + mean   


 


	def evaluate(self): 

		self.fitness[0] = self.fit_func(self.population[0,:])   
		self.best_fit = self.fitness[0] 

		for i in range(self.pop_size):
			self.fitness[i] = self.fit_func(self.population[i,:]) 
			if (self.best_fit> self.fitness[i]):  
				self.best_fit =  self.fitness[i]
				self.best_index = i

			self.num_evals = self.num_evals + 1

		


 
   # calculates the magnitude of a vector
	def mod(self, List):
		sum = 0
		for i in range(self.dimen):
			sum += (List[i] * List[i] )
		return np.sqrt(sum)
 
 

	def parent_centric_xover(self, current):
		centroid = np.zeros(self.dimen)
		tempar1 = np.zeros(self.dimen)
		tempar2 = np.zeros(self.dimen)

		temp_rand = np.zeros(self.dimen)

		d = np.zeros(self.dimen)
		D = np.zeros(self.num_parents)

		temp1, temp2, temp3 = (0,0,0)

		diff = np.zeros((self.num_parents, self.dimen))
		for i in range(self.dimen):
			for u in range(self.num_parents): 
				#print self.temp_index[u], self.population[self.temp_index[u],i], ' self.population[self.temp_index[u],i]'
				centroid[i]  = centroid[i] +  self.population[self.temp_index[u],i]
		centroid   = centroid / self.num_parents 
 
 

  # calculate the distace (d) from centroid to the index parent self.temp_index[0]
  # also distance (diff) between index and other parents are computed

		for j in range(1, self.num_parents):
			for i in range(self.dimen):
				if j == 1:
					d[i]= centroid[i]  - self.population[self.temp_index[0],i]
				diff[j, i] = self.population[self.temp_index[j], i] - self.population[self.temp_index[0],i]


			if (self.mod(diff[j,:]) < self.EPSILON):
				print 'Points are very close to each other. Quitting this run' 
				#return 0 

 



		dist = self.mod(d) 

		#print dist, '                                         ---- dist --+++++++++++++++++++-------- '

		if (dist < self.EPSILON):
			print " Error -  points are very close to each other. Quitting this run   "
			return 0

		# orthogonal directions are computed
		for j in range(1, self.num_parents):
			temp1 = self.inner(diff[j,:] , d )

			if ((self.mod(diff[j,:]) * dist) == 0):
				print "Division by zero"
				temp2 = temp1 / (1)
			else:
				temp2 = temp1 / (self.mod(diff[j,:]) * dist)
 

			temp3 = 1.0 - np.power(temp2, 2)

			D[j] = self.mod(diff[j]) * np.sqrt(np.abs(temp3)) 


 
 
		D_not = 0.0
		for i in range(1, self.num_parents):
			D_not += D[i]
		D_not /= (self.num_parents - 1) # this is the average of the perpendicular distances from all other parents (minus the index parent) to the index vector

		#for j in range(self.dimen):
		#	tempar1[j] =  self.rand_normal(0, D_not * self.sigma_zeta) 

		tempar1 = np.random.normal(0,  self.sigma_eta * D_not , self.dimen) #rand_normal(0, D_not * sigma_eta); 

		tempar2 = tempar1 
 
		if(np.power(dist, 2) == 0):
			print " division by zero: part 2"
			tempar2  = tempar1  
		else:
			tempar2  = tempar1  - (    np.multiply(self.inner(tempar1, d) , d )  ) / np.power(dist, 2.0) 
 
		tempar1= tempar2   
 
		self.sub_pop[current,:] = self.population[self.temp_index[0],:] + tempar1  
		print current, self.sub_pop[current,:], ' * new pcx'


		#for j in range(self.dimen):
		#	temp_rand[j] =  self.rand_normal(0, self.sigma_zeta) 
  
 
		temp_rand = np.random.normal(0, self.sigma_zeta, self.dimen) 


		self.sub_pop[current,:] += np.multiply(temp_rand ,  d ) 
		print current, self.sub_pop[current,:], ' new pcx'
  
  # the child is included in the newpop and is evaluated
		self.sp_fit[current] = self.fit_func(self.sub_pop[current,:])
		self.num_evals += 1 

		return 1


	def inner(self, ind1, ind2):

		sum = 0.0
		for i in range(self.dimen):
			sum += (ind1[i] * ind2[i] )
		return  sum 
 
 
 

	def sort_population(self):

		dbest = 99
 
		for i in range(self.children + self.family):
			self.list[i] = i

		print self.list, ' sort_population list'

		for i in range(self.children + self.family - 1):
			dbest = self.sp_fit[self.list[i]] 



			for j in range(i + 1, self.children + self.family):

				if(self.sp_fit[self.list[j]]  < dbest):
					dbest = self.sp_fit[self.list[j]] 
					temp = self.list[j]
					self.list[j] = self.list[i]
					self.list[i] = temp

 



	def replace_parents(self): #here the best (1 or 2) individuals replace the family of parents

		for j in range(self.family):
			self.population[ self.parents[j],:]  =  self.sub_pop[ self.list[j],:] # Update population with new species 


			fx = self.fit_func(self.population[ self.parents[j],:]) 

			print self.parents[j],  self.population[ self.parents[j],:] , fx,'                               replaced                        @@'
 

			self.fitness[self.parents[j]]   =  fx #self.sp_fit[self.list[j]]
			self.num_evals += 1



	def family_members(self): #//here a random family (1 or 2) of parents is created who would be replaced by good individuals

		swp = 0


		for i in range(self.pop_size): 
			self.parents[i] = i

		for i in range(self.family):
			randomIndex = random.randint(0, self.pop_size - 1) + i # Get random index in population

			if randomIndex > (self.pop_size-1):
				randomIndex = self.pop_size-1

			#print randomIndex, '               randomIndex   --------------------------------  ++++++++++++++   '
			swp = self.parents[randomIndex]
			self.parents[randomIndex] = self.parents[i]
			self.parents[i] = swp 





	def find_parents(self): #here the parents to be replaced are added to the temporary subpopulation to assess their goodness against the new solutions formed which will be the basis of whether they should be kept or not 


		self.family_members()

		for j in range(self.family):
			self.sub_pop[self.children + j, :] = self.population[self.parents[j],:]


			fx = self.fit_func(self.sub_pop[self.children + j, :])     
			self.sp_fit[self.children + j]  = fx


			#print self.sub_pop[self.children + j, :] , self.children + j, fx, j , self.parents[j], '                             								******* -----------@' 
			

			 
			self.num_evals += 1 

 
  
 

	def best_inpopulation(self): 
		self.best_fit = self.fitness[0]
		for y in range(self.pop_size):
			if  self.fitness[y]< best_fit:
				self.best_index = y
				self.best_fit = self.fitness[y]
 

	def worst_inpopulation(self):  
		worst_fit = 0
		for y in range(self.pop_size): 
			if  self.fitness[y] > worst_fit:
				self.worst_index = y
				self.worst_fit = self.fitness[y]
 

	'''def order_population(self): 

		for y in range( self.pop_size - 1):
			min_idx = y
			for x in range(y + 1, self.pop_size):
				if self.fitness[x] < self.fitness[min_idx]:
					min_idx = x
				if (self.rank[y] < self.rank[min_idx]):  # Only swap ranks if next rank is less else leave as is
					tempRank = self.rank[min_idx]   # Take previous best index rank
					self.rank[min_idx]  = self.rank[y]   # Swap ranks
					self.rank[y]  = tempRank
 
		for x in range(self.pop_size):
			if self.rank[x] == 0:
				self.best_index = x
				self.best_fit = self.fitness[x]'''

 

	def random_parents(self ):
 
		for i in range(self.pop_size): 
			self.temp_index[i] = i

		swp=self.temp_index[0] 
		self.temp_index[0]=self.temp_index[self.best_index] 

		self.temp_index[self.best_index]  = swp 

		 #best is always included as a parent and is the index parent
		  # this can be changed for solving a generic problem

		for i in range(1, self.rand_parents): 
			index= np.random.randint(self.pop_size)+i  

			if index > (self.pop_size-1):
				index = self.pop_size-1
			swp=self.temp_index[index] 
			self.temp_index[index]=self.temp_index[i] 
			self.temp_index[i]=swp 


 
 



	def evolve(self, outfile   ):

		#np.savetxt(outfile, self.population, fmt = '%1.2f' )

		pop = np.loadtxt("pop.txt" )
		genIndex = np.loadtxt("out3.txt" )
		mom = np.loadtxt("out2.txt" )

		self.population = pop 

		#print genIndex, ' genIndex'
		print self.population

		 
 
		tempfit = 0

		prevfitness = 99  

		self.evaluate() 

		print self.fitness

		tempfit= self.fitness[self.best_index]

		print tempfit, self.best_fit, self.best_index, ' initial best fit  -------------------------------------------------------------------'



		vxx = 0



		print self.best_fit, ' is initial best fit   ------------ ++ ------------ '


		while(self.num_evals < self.max_evals):  


			#print self.sub_pop, ' --------------------- initial sub_pop'
			#print self.sp_fit, ' --------------------- initial sp_fit'

			tempfit = self.best_fit

			self.random_parents()


			print self.temp_index  , '  -------------------- --------------------------- index of rand_parents'

			#print self.temp_index, ' self.temp'

			#print vxx, '  i ****                            x'  

			#print self.temp_index, ' temp_index'


			for i in range(self.children):
				#print i , '    *                        i     -------------->>>>>> ----------+++++++++++++++++++++++++++'
				tag = self.parent_centric_xover(i) 
				if (tag == 0):
					break
			#if (tag == 0):
				#	break




			self.find_parents()


			self.parents = mom[vxx,:].astype(int)
 





			vxx= vxx + 1



			print self.parents,  '   -------------------------  ------------------------------------------------ self.parents'

			print self.sub_pop, '  -----------------------find parents sub_pop'
			print self.sp_fit, '------------------------- find parents sp_fit'


			self.sort_population()

			print self.sub_pop, ' sort_population sub_pop'
			print self.sp_fit, ' sort_population sp_fit'




			self.replace_parents()

 


			self.best_index = 0
			tempfit = self.fitness[0]

			print self.sub_pop, '  replace_parents sub_pop'
			print self.sp_fit, ' replace_parents sp_fit'

			 


			for x in range(1, self.pop_size):
				if(self.fitness[x] < tempfit):

					self.best_index = x
					tempfit  =  self.fitness[x]

					#print self.best_fit, x, ' is initial best fit   ------------ ++ ------------ '


			print self.fitness, ' fitness'

			print self.sp_fit, ' sp_fit'


			#print self.num_evals, self.best_index, self.best_fit, ' best so far --------------------  '
			print self.fitness, ' fitness -------------- ** '
			print self.num_evals, self.fitness[self.best_index], self.best_index, ' best sol         ---------------   ** ---------------** --'

			print self.population

			#np.savetxt(outfile, [ self.num_evals, self.best_index, self.best_fit], fmt='%1.5f', newline="\n")


		print self.sub_pop, '  sub_pop'
		print self.population[self.best_index], ' best sol                                         '
		print self.fitness, ' fitness'
 
 
 

def main():


	
	outfile=open('pop_.txt','w')


	 
 
 

	MinCriteria = 0.005  # stop when RMSE reaches MinCriteria ( problem dependent)

	random.seed(time.time())

	max_evals = 200   # need to decide yourself 80000

	pop_size =  20
	num_varibles = 2

	max_limits = np.repeat(5, num_varibles) 
	min_limits = np.repeat(-5, num_varibles) 
 

	g3pcx  = evolution(pop_size, num_varibles, max_evals,  max_limits, min_limits)


	g3pcx.evolve(outfile)
 



 

if __name__ == "__main__": main()