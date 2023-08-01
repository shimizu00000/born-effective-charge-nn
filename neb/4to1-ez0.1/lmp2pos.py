#import numpy as np
import ase.io.vasp
#import pymatgen as mg
#from pymatgen import Lattice, Structure, Molecule
#from pymatgen.io.vasp import Poscar
#from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#from pymatgen.analysis.defects.core import Defect, Vacancy, Interstitial
#from pymatgen.core.sites import PeriodicSite
#from pymatgen.analysis.defects.generators import VacancyGenerator
#from ase.calculators.lammpsrun import write_lammps_data
import ase.io.lammpsdata
import ase.io.lammpsrun

#poscar = Poscar.from_file("POSCAR")
#structure = poscar.structure
#sp = structure.species

#ASE
#atoms = ase.io.vasp.read_vasp("POSCAR")
species = ['Li','P','O']
#write_lammps_data('data.lammps',atoms,specorder=species)

atoms = ase.io.lammpsdata.read_lammps_data("data.lammps",None,"atomic",False,"metal")
ase.io.vasp.write_vasp("POSCAR",atoms,label='POSCAR by ASE',direct=True,sort=True)


#cell = ase.io.vasp.read_vasp("POSCAR331")
#ase.io.vasp.write_vasp("POSCAR331",cell*(3,3,1),label='331supercell',direct=True,sort=True)
#ase.io.vasp.write_vasp("POSCAR331",cell*(3,3,1),label='331supercell',direct=False,sort=True)
#ase.io.vasp.write_vasp("POSCAR331cart",cell*(1,1,1),label='331supercell',direct=False,sort=True)
#ase.io.vasp.write_vasp("POSCAR111",cell*(1,1,1),label='111supercell',direct=False,sort=True)
