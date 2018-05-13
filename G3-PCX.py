# !/usr/bin/python

 
# Rohitash Chandra, Centre for Translational Data Science
# University of Sydey, Sydney NSW, Australia.  2017 c.rohitash@gmail.conm
# https://www.researchgate.net/profile/Rohitash_Chandra


import matplotlib.pyplot as plt
import numpy as np
import random
import time
import math
import os
import shutil




 
class evolution:
	def __init__(self, pop_size, dimen, max_evals,  max_limits, min_limits):
		
		self.EPSILON = 1e-50  # convergence
		self.sigma_eta = 0.1
		self.sigma_zeta = 0.1

		self.children = 2
		self.num_parents = 3 
		self.family = 2
		self.sp_size = self.children + self.family 



		self.population =   np.random.randn( pop_size  , dimen)  * 5  #[SpeciesPopulation(dimen) for count in xrange(pop_size)]
		self.sub_pop =  np.random.randn( self.sp_size , dimen )  * 5  #[SpeciesPopulation(dimen) for count in xrange(NPSize)]

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




	def fit_func(self, x):    # rosenbrock function
		fit = 0
		for j in range(x.size -1): 
			fit += 100.0*(x[j]*x[j] - x[j+1])*(x[j]*x[j] - x[j+1]) + (x[j]-1.0)*(x[j]-1.0)

		return fit

	def evaluate(self): 

		self.fitness[0] = self.fit_func(self.population[0,:])   
		self.best_fit = self.fitness[0] 

		for i in range(self.pop_size):
			self.fitness[i] = self.fit_func(self.population[i,:]) 
			if (self.best_fit> self.fitness[i]):  
				self.best_fit =  self.fitness[i]
				self.best_index = i

			self.num_evals = self.num_evals + 1

		


 
  
	'''def RandomParents(self):
		swp = self.temp_index[0]
		self.temp_index[0] = self.temp_index[self.best_index]
		self.temp_index[self.best_index] = swp

		for i in range(1, self.num_parents):
			index = random.randint(1, self.pop_size - 1)
			swp = self.temp_index[index]
			self.temp_index[index] = self.temp_index[i]
			self.temp_index[i] = swp'''

	def mod(self, List):
		sum = 0
		for i in range(self.dimen):
			sum += (List[i] * List[i] )
		return np.sqrt(sum)
 
 

	def parent_centric_xover(self, current):
		centroid = np.zeros(self.dimen)
		tempar1 = np.zeros(self.dimen)
		tempar2 = np.zeros(self.dimen)
		d = np.zeros(self.dimen)
		D = np.zeros(self.num_parents)

		temp1, temp2, temp3 = (0,0,0)

		diff = np.zeros((self.num_parents, self.dimen))

		for i in range(self.dimen):
			for u in range(self.num_parents):
				centroid[i] += self.population[self.temp_index[u], i]
			centroid[i] /= self.num_parents

		print centroid, ' Centroid'


		for j in range(1, self.num_parents):
			for i in range(self.dimen):
				if j == 1:
					d[i] = centroid[i] - self.population[self.temp_index[0],i]

				if(np.isnan(self.population[self.temp_index[j],i] - self.population[self.temp_index[0] ,i])):
					print 'diff nan'
					diff[j][i] = 1
					return 0
				else:
					diff[j][i] = self.population[self.temp_index[j],i] - self.population[self.temp_index[0],i]

			if (self.mod(diff[j])  < self.EPSILON):
				# print "RUN Points are very close to each other. Quitting this run   "
				return 0

		#print d, ' d '


		dist = self.mod(d)

		#print dist, ' dist'

		if (dist < self.EPSILON):
			# print "RUN Points are very close to each other. Quitting this run   "
			return 0

		# orthogonal directions are computed
		for j in range(1, self.num_parents):
			temp1 = np.inner(diff[j].tolist(), d.tolist())

			if ((self.mod(diff[j]) * dist) == 0):
				print "Division by zero"
				temp2 = temp1 / (1)
			else:
				temp2 = temp1 / (self.mod(diff[j]) * dist)
 

			temp3 = 1.0 - np.power(temp2, 2)

			D[j] = self.mod(diff[j]) * np.sqrt(np.abs(temp3))

		D_not = 0
		for i in range(1, self.num_parents):
			D_not += D[i]
		D_not /= (self.num_parents - 1) # this is the average of the perpendicular distances from all other parents (minus the index parent) to the index vector

		#for j in range(self.dimen):
		tempar1 = np.random.normal(0, self.sigma_eta, self.dimen) #rand_normal(0, D_not * sigma_eta);
		tempar2 = tempar1 

		for j in range(self.dimen):
			if(np.power(dist, 2) == 0):
				print " division by zero: part 2"
				tempar2[j] = tempar1[j] - ((np.inner(tempar1, d) * d[j]) / 1)
			else:
				tempar2[j] = tempar1[j] - ((np.inner(tempar1, d) * d[j]) / np.power(dist, 2.0))

		#for j in range(self.dimen):
		tempar1= tempar2 

		#for k in range(self.dimen):
		self.sub_pop[current,:] = self.population[self.temp_index[0],:] + tempar1 

		temp_rand = np.random.normal(0, self.sigma_zeta, 1)

		for k in range(self.dimen):
			self.sub_pop[current,k] += temp_rand  #* d[k])
 

 
 



		return 1

	def family_members(self):

		swp = 0

		for i in range(self.family):
			randomIndex = random.randint(0, self.pop_size - 1)  # Get random index in population
			swp = self.parents[randomIndex]
			self.parents[randomIndex] = self.parents[i]
			self.parents[i] = swp

  
 

	def sort_population(self):

		dbest = 99
 
		for i in range(self.children + self.family):
			self.list[i] = i

		for i in range(self.children + self.family - 1):
			dbest = self.sp_fit[self.list[i]] 

			for j in range(i + 1, self.children + self.family):

				if(self.sp_fit[self.list[j]]  < dbest):
					dbest = self.sp_fit[self.list[j]] 
					temp = self.list[j]
					self.list[j] = self.list[i]
					self.list[i] = temp


	def replace_parents(self):

		for j in range(self.family):
			self.population[ self.parents[j],:]  =  self.sub_pop[ self.list[j],:] # Update population with new species 
 
			fitness = self.fit_func(self.population[j,:]) 
			self.fitness[self.parents[j]]   =  fitness #self.sp_fit[self.list[j]]
			self.num_evals += 1

	def find_parents(self):

		self.family_members()

		for j in range(self.family):
			self.sub_pop[self.children + j, :] = self.population[self.parents[j],:]
			fitness = self.fit_func(self.sub_pop[self.children + j, :])     
			self.sp_fit[self.children + j]  = fitness

			print j, fitness, ' find_parents fit'
			self.num_evals += 1
			#self.sp_fit[self.children + j] = self.fitness[self.parents[j]]

 

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
 

	def order_population(self): 

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
				self.best_fit = self.fitness[x]


		print('Completed Population Sort by Fitness') 


	def random_parents(self):
 
		for i in range(self.pop_size): 
			self.temp_index[i] = i

		swp=self.temp_index[0] 
		self.temp_index[0]=self.temp_index[self.best_index] 

		for i in range(self.rand_parents): 
			index=  np.random.randint(self.pop_size)+i 
			print index, ' index ..'

			if index > (self.pop_size-1):
				index = self.pop_size-1
				swp=self.temp_index[index] 
				self.temp_index[index]=self.temp_index[i] 
				self.temp_index[i]=swp 



 



	def evolve(self   ):
 
		tempfit = 99

		prevfitness = 99  

		self.evaluate() 

		tempfit= self.fitness[self.best_index]


		while(self.num_evals < self.max_evals):  

			tempfit = self.best_fit

			self.random_parents()

			#print self.temp_index, ' temp_index'

			for i in range(self.children):
				tag = self.parent_centric_xover(i) 
				if (tag == 0):
					break
			if (tag == 0):
					break
 


			self.find_parents()
			self.sort_population()
			self.replace_parents()

			
			print self.sp_fit, ' sp_fit'
 

			self.best_index = 0
			self.best_fit = self.fitness[0]


			for i in range(1, self.pop_size):
				if(self.fitness[i] < self.best_fit):

					self.best_index = i
					self.best_fit = tempfit

			print self.num_evals, self.best_index, self.best_fit, ' best so far'

		print self.sub_pop, '  sub_pop'
		print self.population[self.best_index], ' best sol'
		print self.fitness, ' fitness'
 
 
 

def main():

	 
 
 

	MinCriteria = 0.005  # stop when RMSE reaches MinCriteria ( problem dependent)

	random.seed(time.time())

	max_evals = 2000  # need to decide yourself 80000

	pop_size = 100
	num_varibles = 5

	max_limits = [1, 1, 1, 1, 1]
	min_limits = [0, 0, 0, 0, 0]
 

	g3pcx  = evolution(pop_size, num_varibles, max_evals,  max_limits, min_limits)


	g3pcx.evolve()
 



 

if __name__ == "__main__": main()
