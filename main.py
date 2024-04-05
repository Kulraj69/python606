import random
import sys
import math
import copy
import os

# Output directory
output_dir = "./output"

# Simulation Parameters
N = 500
sigma = 1
epsilon = 1
trunc = 3 * sigma
truncsq = trunc ** 2
steps = 100000
temp = 8.5e-1
density = 1.0e-3
L = (N / density) ** (1.0 / 3.0)
halfL = L / 2
particles = []

# Some helper functions
def wrap(particle):
    '''Apply periodic boundary conditions.'''
    for i in range(3):
        if particle[i] > L:
            particle[i] -= L
        elif particle[i] < 0:
            particle[i] += L
    return particle

def distancesq(particle1, particle2):
    '''Gets the squared distance between two particles, applying the minimum image convention.'''
    dx = particle1[0] - particle2[0]
    dy = particle1[1] - particle2[1]
    dz = particle1[2] - particle2[2]
    # Minimum image convention
    dx -= L * round(dx / L)
    dy -= L * round(dy / L)
    dz -= L * round(dz / L)
    return dx ** 2 + dy ** 2 + dz ** 2

def energy(particles):
    '''Gets the energy of the system'''
    energy = 0
    for particle1 in range(0, len(particles) - 1):
        for particle2 in range(particle1 + 1, len(particles)):
            dist = distancesq(particles[particle1], particles[particle2])
            if dist <= truncsq:
                energy += 4 * (1 / dist ** 6) - (1 / dist ** 3)
    return energy

def particleEnergy(particle, particles, p):
    '''Gets the energy of a single particle.'''
    part_energy = 0
    for i, particle2 in enumerate(particles):
        if i != p:
            dist = distancesq(particle, particle2)
            if dist <= truncsq:
                part_energy += 4 * (1 / dist ** 6) - (1 / dist ** 3)
    return part_energy

def writeEnergy(step, en):
    '''Writes the energy to a file.'''
    with open('./output/energy.csv', 'a') as f:
        f.write('{0},{1}\n'.format(step, en))

if __name__ == "__main__":
    # Checking output directory exists
    if not os.path.exists(output_dir):
        # If it doesn't exist, create the directory
        os.makedirs(output_dir)

    # Clear files if they already exist.
    if os.path.exists('./output/energy.csv'):
        os.remove('./output/energy.csv')

    # Initiating csv file:
    with open('./output/energy.csv', 'a') as f:
        f.write('microstate,energy\n')

    # Initialize the simulation box:
    for particle in range(0, N):
        x_coord = random.uniform(0, L)
        y_coord = random.uniform(0, L)
        z_coord = random.uniform(0, L)
        particles.append([x_coord, y_coord, z_coord])

    # Calculate initial energy
    en = energy(particles)

    # MC
    for step in range(0, steps):
        sys.stdout.write("\rStep: {0} Energy: {1}".format(step, en))
        sys.stdout.flush()

        # Choose a particle to move at random.
        p = random.randint(0, N - 1)

        # Move particle and evaluate energy
        this_particle = copy.deepcopy(particles[p])
        prev_E = particleEnergy(this_particle, particles, p)
        for i in range(3):
            this_particle[i] += random.uniform(-1, 1)
        this_particle = wrap(this_particle)
        new_E = particleEnergy(this_particle, particles, p)
        deltaE = new_E - prev_E

        # Acceptance rule enforcing Boltzmann statistics
        if deltaE < 0:
            particles[p] = this_particle
            en += deltaE
        else:
            rand = random.random()
            if math.exp(-deltaE / temp) > rand:
                particles[p] = this_particle
                en += deltaE
        writeEnergy(str(step), str(en))
