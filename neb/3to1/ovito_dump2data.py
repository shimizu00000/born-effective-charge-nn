#
# conda activate ovito3
#
import ovito
from ovito.io import import_file, export_file

pipeline = import_file("dump.melt")
pipeline.add_to_scene()

def setup_particle_types(frame,data):
    types = data.particles_.particle_types_
    types.type_by_id_(1).mass = 7.010
    types.type_by_id_(2).mass = 30.974
    types.type_by_id_(3).mass = 32.066
#1 7.010 Li
#2 30.974 P
#3 32.066 S

pipeline.modifiers.append(setup_particle_types)

export_file(pipeline, "data_all.lammps", "lammps/data", multiple_frames=True)
export_file(pipeline, "data.*.lammps", "lammps/data", multiple_frames=True)
