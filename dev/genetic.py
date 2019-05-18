from __future__ import division
import os
import glob
import string
import random
import numpy as np
from scipy.odr import *
import json
import matplotlib.pyplot as plt
import CoolProp
import CoolProp.CoolProp as CP

LIBRARY = [i / 6.0 for i in range(1, 151)] + [0.35 + i / 2000 for i in range(1, 100)] + [0.05 + 0.001 * i for i in range(1, 100)] + [i + 0.5 for i in range(10)]
#LIBRARY = [i/1000 for i in range(1,20000)]


class Sample(object):
    def __init__(self, v):
        self.v = v


class GeneticAncillaryFitter(object):
    def __init__(self,
               num_samples=600,  # Have this many chromos in the sample group
               num_selected=60,  # Have this many chromos in the selected group
               mutation_factor=2,  # Randomly mutate 1/n of the chromosomes
               num_powers=5,  # How many powers in the fit
               Ref='R407C',
               value='rhoV',
               addTr=True,
               values=None,
               Tlims=None
                ):
        self.num_samples = num_samples
        self.num_selected = num_selected
        self.mutation_factor = mutation_factor
        self.num_powers = num_powers
        self.addTr = addTr
        self.value = value
        self.Ref = Ref

        # Thermodynamics
        from CoolProp.CoolProp import PropsSI

        if values is None:
            self.Tc = PropsSI(Ref, 'Tcrit')
            self.pc = PropsSI(Ref, 'pcrit')
            self.rhoc = PropsSI(Ref, 'rhomolar_critical')
            self.Tmin = PropsSI(Ref, 'Tmin')
            if Tlims is None:
                self.T = np.append(np.linspace(self.Tmin + 1e-14, self.Tc - 1, 150), np.logspace(np.log10(self.Tc - 1), np.log10(self.Tc) - 1e-15, 40))
            else:
                self.T = np.linspace(Tlims[0], Tlims[1])
            self.pL = np.array(PropsSI('P', 'T', self.T, 'Q', [0] * len(self.T), Ref))
            self.pV = np.array(PropsSI('P', 'T', self.T, 'Q', [1] * len(self.T), Ref))
            self.rhoL = PropsSI('Dmolar', 'T', self.T, 'Q', [0] * len(self.T), Ref)
            self.rhoV = PropsSI('Dmolar', 'T', self.T, 'Q', [1] * len(self.T), Ref)
        else:
            self.Tc = values['Tcrit']
            self.pc = values['pcrit']
            self.rhoc = values['rhocrit']
            self.Tmin = values['Tmin']
            self.T = np.array(values['T'])
            self.p = np.array(values['p'])
            self.pL = np.array(values['p'])
            self.pV = np.array(values['p'])
            self.rhoL = np.array(values['rhoL'])
            self.rhoV = np.array(values['rhoV'])

        self.logpLpc = (np.log(self.pL) - np.log(self.pc))
        self.logpVpc = (np.log(self.pV) - np.log(self.pc))
        self.rhoLrhoc = np.array(self.rhoL) / self.rhoc
        self.rhoVrhoc = np.array(self.rhoV) / self.rhoc
        self.logrhoLrhoc = np.log(self.rhoL) - np.log(self.rhoc)
        self.logrhoVrhoc = np.log(self.rhoV) - np.log(self.rhoc)

        self.x = 1.0 - self.T / self.Tc

        MM = PropsSI(Ref, 'molemass')
        self.T_r = self.Tc

        if self.value == 'pL':
            self.LHS = self.logpLpc.copy()
            self.EOS_value = self.pL.copy()
            if self.addTr == False:
                self.description = "p' = pc*exp(sum(n_i*theta^t_i))"
            else:
                self.description = "p' = pc*exp(Tc/T*sum(n_i*theta^t_i))"
            self.reducing_value = self.pc
        elif self.value == 'pV':
            self.LHS = self.logpVpc.copy()
            self.EOS_value = self.pV.copy()
            if self.addTr == False:
                self.description = "p'' = pc*exp(sum(n_i*theta^t_i))"
            else:
                self.description = "p'' = pc*exp(Tc/T*sum(n_i*theta^t_i))"
            self.reducing_value = self.pc
        elif self.value == 'rhoL':
            self.LHS = self.logrhoLrhoc.copy()
            self.EOS_value = self.rhoL
            if self.addTr == False:
                self.description = "rho' = rhoc*exp(sum(n_i*theta^t_i))"
            else:
                self.description = "rho' = rhoc*exp(Tc/T*sum(n_i*theta^t_i))"
            self.reducing_value = self.rhoc
        elif self.value == 'rhoV':
            self.LHS = self.logrhoVrhoc.copy()
            self.EOS_value = self.rhoV
            if self.addTr == False:
                self.description = "rho'' = rhoc*exp(sum(n_i*theta^t_i))"
            else:
                self.description = "rho'' = rhoc*exp(Tc/T*sum(n_i*theta^t_i))"
            self.reducing_value = self.rhoc
        elif self.value == 'rhoLnoexp':
            self.LHS = (self.rhoLrhoc - 1).copy()
            self.EOS_value = self.rhoL
            self.description = "rho' = rhoc*(1+sum(n_i*theta^t_i))"
            self.reducing_value = self.rhoc
        else:
            raise ValueError

        if self.value == 'rhoLnoexp' and self.addTr:
            raise ValueError('Invalid combination')

        if self.addTr:
            self.LHS *= self.T / self.Tc

    def generate_random_chromosomes(self,):
        '''
        Create a list of random chromosomes to seed our algorithm.
        '''
        chromos = []
        while len(chromos) < self.num_samples:
            chromos.append(Sample(sorted(random.sample(LIBRARY, self.num_powers))))
        return chromos

    def fitness(self, chromo):
        '''
        Fitness of a chromo is the sum of the squares of the error of the correlation
        '''

        # theta^t where the i-th row is equal to theta^t[i]
        # We need these terms later on to build the A and b matrices
        theta_t = (self.x.reshape(-1, 1)**chromo.v).T

        # TODO: more numpy broadcasting should speed this up even more
        # Another few times faster ought to be possible
        I = len(chromo.v)
        A = np.zeros((I, I))
        b = np.zeros((I, 1))
        for i in range(I):
            for j in range(I):
                A[i, j] = np.sum(theta_t[i] * theta_t[j])
            b[i] = np.sum(theta_t[i] * self.LHS)

        # If you end up with a singular matrix, quit this run
        try:
            n = np.linalg.solve(A, b).T
        except np.linalg.linalg.LinAlgError as E:
            chromo.fitness = 1e99
            return

        chromo.beta = n
        RHS = np.sum(n * self.x.reshape(-1, 1)**chromo.v, axis=1)

        if self.addTr:
            RHS *= self.Tc / self.T

        if self.value in ['pL', 'pV']:
            fit_value = np.exp(RHS) * self.pc
        elif self.value in ['rhoL', 'rhoV']:
            fit_value = np.exp(RHS) * self.rhoc
        elif self.value == 'rhoLnoexp':
            fit_value = self.rhoc * (1 + RHS)
        else:
            raise ValueError

        max_abserror = np.max(np.abs((fit_value / self.EOS_value) - 1) * 100)

        chromo.fitness = max_abserror
        chromo.fit_value = fit_value
        chromo.max_abserror = max_abserror

        return chromo.fitness

    def tourny_select_chromo(self, samples):
        '''
        Randomly select two chromosomes from the samples, then return the one
        with the best fitness.
        '''
        a = random.choice(samples)
        b = random.choice(samples)
        if a.fitness < b.fitness:
            return a
        else:
            return b

    def breed(self, a, b):
        '''
        Breed two chromosomes by splicing them in a random spot and combining
        them together to form two new chromos.
        '''
        splice_pos = random.randrange(len(a.v))
        new_a = a.v[:splice_pos] + b.v[splice_pos:]
        new_b = b.v[:splice_pos] + a.v[splice_pos:]
        return Sample(sorted(new_a)), Sample(sorted(new_b))

    def mutate(self, chromo):
        '''
        Mutate a chromosome by changing one of the parameters, but only if it improves the fitness
        '''
        v = chromo.v
        if hasattr(chromo, 'fitness'):
            old_fitness = chromo.fitness
        else:
            old_fitness = self.fitness(chromo)

        for i in range(10):
            pos = random.randrange(len(chromo.v))
            chromo.v[pos] = random.choice(LIBRARY)
            new_fitness = self.fitness(chromo)
            if new_fitness < old_fitness:
                return chromo
            else:
                return Sample(sorted(v))

    def run(self):
        # Create a random sample of chromos
        samples = self.generate_random_chromosomes()

        # Calculate the fitness for the initial chromosomes
        for chromo in samples:
            self.fitness(chromo)
#         print '#'

        decorated = sorted([(sample.fitness, sample) for sample in samples])
        samples = [s for sv, s in decorated]
        values = [sv for sv, s in decorated]
        plt.plot(values[0:len(values) // 2])
        plt.close()

        # Main loop: each generation select a subset of the sample and breed from
        # them.
        generation = -1
        while generation < 0 or samples[0].fitness > 0.02 or (generation < 3 and generation < 15):
            generation += 1

            # Generate the selected group from sample- take the top 10% of samples
            # and tourny select to generate the rest of selected.
            ten_percent = int(len(samples) * .1)
            selected = samples[:ten_percent]
            while len(selected) < self.num_selected:
                selected.append(self.tourny_select_chromo(samples))

            # Generate the solution group by breeding random chromos from selected
            solution = []
            while len(solution) < self.num_samples:
                solution.extend(self.breed(random.choice(selected),
                                         random.choice(selected)))

            # Apply a mutation to a subset of the solution set
            mutate_indices = random.sample(range(len(solution)), len(solution) // self.mutation_factor)
            for i in mutate_indices:
                solution[i] = self.mutate(solution[i])

            for chromo in solution:
                self.fitness(chromo)
#             print '#'

            decorated = sorted([(sample.fitness, sample) for sample in solution])
            samples = [s for sv, s in decorated]

#             print '------------------  Top 10 values  ---------------'
#             for sample in samples[0:10]:
#                 print sample.v, sample.fitness, sample.max_abserror

#             print '// Max error is ',samples[0].max_abserror,'% between',np.min(self.T),'and',np.max(self.T),'K'
#             print str(samples[0].v), samples[0].beta.tolist()

            # Print useful stats about this generation
            (min, median, max) = [samples[0].fitness, samples[len(samples) // 2].fitness, samples[-1].fitness]
#             print("{0} best value: {1}. fitness: best {2}, median {3}, worst {4}".format(generation, samples[0].v, min, median, max))

            # If the string representations of all the chromosomes are the same, stop
            if len(set([str(s.v) for s in samples[0:5]])) == 1:
                break

        self.fitness(samples[0])

        print(self.value)
        print('// Max error is ' + str(samples[0].max_abserror) + ' % between ' + str(np.min(self.T)) + ' and ' + str(np.max(self.T)) + ' K')

        self.fit_value = samples[0].fit_value

        j = dict()
        j['n'] = samples[0].beta.squeeze().tolist()
        j['t'] = samples[0].v
        j['Tmin'] = np.min(self.T)
        j['Tmax'] = np.max(self.T)
        j['type'] = self.value
        j['using_tau_r'] = self.addTr
        j['reducing_value'] = self.reducing_value
        j['T_r'] = self.Tc

        # Informational, not used
        j['max_abserror_percentage'] = samples[0].max_abserror
        j['description'] = self.description

        return j


def build_ancillaries(name, **kwargs):

    j = dict()
    j['ANCILLARIES'] = dict()
    gaf = GeneticAncillaryFitter(Ref=name, value='pL', addTr=True, num_powers=6, **kwargs)
    j['ANCILLARIES']['pL'] = gaf.run()
    gaf = GeneticAncillaryFitter(Ref=name, value='pV', addTr=True, num_powers=6, **kwargs)
    j['ANCILLARIES']['pV'] = gaf.run()
    gaf = GeneticAncillaryFitter(Ref=name, value='rhoLnoexp', addTr=False, num_powers=6, **kwargs)
    j['ANCILLARIES']['rhoL'] = gaf.run()
    gaf = GeneticAncillaryFitter(Ref=name, value='rhoV', addTr=True, num_powers=6, **kwargs)
    j['ANCILLARIES']['rhoV'] = gaf.run()

    fp = open(os.path.join('ancillaries', name + '_anc.json'), 'w')
    print(json.dumps(j, indent=2), file=fp)
    fp.close()


def build_all_ancillaries():
    for fluid in sorted(CoolProp.__fluids__):
        print(fluid)
        if fluid in ['SES36']:
            build_ancillaries(fluid, Tlims=[CP.Props(fluid, 'Ttriple'), CP.Props(fluid, 'Tcrit') - 1])
        elif fluid == 'R507A':
            build_ancillaries(fluid, Tlims=[CP.Props(fluid, 'Ttriple'), CP.Props(fluid, 'Tcrit') - 0.1])
        elif fluid == 'R407F':
            build_ancillaries(fluid, Tlims=[CP.Props(fluid, 'Ttriple'), CP.Props(fluid, 'Tcrit') - 2])
        else:
            build_ancillaries(fluid)


if __name__ == "__main__":

    fluid = 'Methanol'
    RPfluid = fluid
    build_ancillaries(RPfluid, Tlims=[CP.PropsSI(fluid, 'Ttriple'), CP.PropsSI(fluid, 'Tcrit') - 0.01])

    #~ build_all_ancillaries()
#     inject_ancillaries()
