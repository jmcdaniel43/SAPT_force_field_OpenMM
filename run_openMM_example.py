from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
#******** this is module that goes with sapt force field files to generate exclusions
from sapt_exclusions import *
#***************************

temperature=300*kelvin
pressure = 1.0*atmosphere
barofreq = 100
pdb = PDBFile('1ions_400water.pdb')
strdir = 'simulation_output/'

integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
integ_md.setMaxDrudeDistance(0.02)  # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)

pdb.topology.loadBondDefinitions('sapt_residues.xml')
pdb.topology.createStandardBonds();
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('sapt.xml')
modeller.addExtraParticles(forcefield)

system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.1*nanometer, constraints=None, rigidWater=True)
nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
nbondedForce.setNonbondedMethod(NonbondedForce.PME)
customNonbondedForce.setNonbondedMethod(min(nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))
customNonbondedForce.setUseLongRangeCorrection(True)

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

barostat = MonteCarloBarostat(pressure,temperature,barofreq)
system.addForce(barostat)
barofreq = barostat.getFrequency()

platform = Platform.getPlatformByName('CPU')
#platform = Platform.getPlatformByName('OpenCL')
#properties = {'OpenCLPrecision': 'mixed'}

simmd = Simulation(modeller.topology, system, integ_md, platform)
#simmd = Simulation(modeller.topology, system, integ_md, platform, properties)
simmd.context.setPositions(modeller.positions)


#************************************************
#         IMPORTANT: generate exclusions for SAPT-FF
#
sapt_exclusions = sapt_generate_exclusions(simmd,system,modeller.positions)
#
#************************************************

# initial energies
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))
for j in range(system.getNumForces()):
    f = system.getForce(j)
    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))


simmd.reporters = []
simmd.reporters.append(DCDReporter(strdir+'md_npt.dcd', 1000))
simmd.reporters.append(CheckpointReporter(strdir+'md_npt.chk', 10000))
simmd.reporters[1].report(simmd,state)

print('Simulating...')

for i in range(1,10):
    simmd.step(5)
    print(i,strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    print(i,datetime.now())
    state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))
    for j in range(system.getNumForces()):
        f = system.getForce(j)
        print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))

print('Done!')

exit()
