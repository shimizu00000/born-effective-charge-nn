import ase.io.vasp
import ase.io.lammpsdata
import ase.io.lammpsrun

#ASE
species = ['Li','P','O']
atoms = ase.io.lammpsdata.read_lammps_data("data.lammps",None,"atomic",False,"metal")
ase.io.vasp.write_vasp("POSCAR",atoms,label='POSCAR by ASE',direct=True,sort=True)


#cell = ase.io.vasp.read_vasp("POSCAR331")
#ase.io.vasp.write_vasp("POSCAR331",cell*(3,3,1),label='331supercell',direct=True,sort=True)
#ase.io.vasp.write_vasp("POSCAR331",cell*(3,3,1),label='331supercell',direct=False,sort=True)
#ase.io.vasp.write_vasp("POSCAR331cart",cell*(1,1,1),label='331supercell',direct=False,sort=True)
#ase.io.vasp.write_vasp("POSCAR111",cell*(1,1,1),label='111supercell',direct=False,sort=True)
