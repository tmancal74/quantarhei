# -*- coding: utf-8 -*-

import random
import numpy

from deap import base
from deap import creator
from deap import tools, algorithms

creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
# Attribute generator
toolbox.register("attr_float", random.uniform, 1.0, 5.0)
# Structure initializers
toolbox.register("individual", tools.initRepeat, creator.Individual, 
    toolbox.attr_float, 3)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def evalOneMax(ind):
    dt = 0.1
    tt = numpy.array([x*dt for x in range(0,1000)])
    a = 2.0
    b = 50.0
    c = 5.0
    yy = a*numpy.exp(-((tt-b)/c)**2)
    zz = ind[0]*numpy.exp(-((tt-ind[1])/ind[2])**2)
    ret = sum(numpy.abs(yy - zz))/(len(yy)*dt)
    print(ret)
    return ret,
    
    
toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)


    
def main():
    pop = toolbox.population(n=30)
    hof = tools.HallOfFame(1, similar=numpy.allclose)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)
    
    pop, log = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.5, ngen=40, 
                                   stats=stats, halloffame=hof, verbose=True)    
#def main():   
#    pop = toolbox.population(n=10)
#    
#    # Evaluate the entire population
#    fitnesses = list(map(toolbox.evaluate, pop))
#    for ind, fit in zip(pop, fitnesses):
#        ind.fitness.values = fit
#
#    NGEN = 20
#    CXPB = 0.5
#    MUTPB = 0.5
#    
#    # Begin the evolution
#    for g in range(NGEN):
#        print("-- Generation %i --" % g)
#
#        # Select the next generation individuals
#        offspring = toolbox.select(pop, len(pop))
#        # Clone the selected individuals
#        offspring = list(map(toolbox.clone, offspring))
#        
#        # Apply crossover and mutation on the offspring
#        for child1, child2 in zip(offspring[::2], offspring[1::2]):
#            if random.random() < CXPB:
#                toolbox.mate(child1, child2)
#                del child1.fitness.values
#                del child2.fitness.values
#
#        for mutant in offspring:
#            if random.random() < MUTPB:
#                toolbox.mutate(mutant)
#                del mutant.fitness.values
#
#        # Evaluate the individuals with an invalid fitness
#        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
#        fitnesses = map(toolbox.evaluate, invalid_ind)
#        for ind, fit in zip(invalid_ind, fitnesses):
#            ind.fitness.values = fit
#            
#        pop[:] = offspring
#        
#        # Gather all the fitnesses in one list and print the stats
#        fits = [ind.fitness.values[0] for ind in pop]
#        
#        length = len(pop)
#        mean = sum(fits) / length
#        sum2 = sum(x*x for x in fits)
#        std = abs(sum2 / length - mean**2)**0.5
#        
#        print("  Min %s" % min(fits))
#        print("  Max %s" % max(fits))
#        print("  Avg %s" % mean)
#        print("  Std %s" % std)
#        
    print(pop)
        
if __name__ == "__main__":
    
    main()
    
    