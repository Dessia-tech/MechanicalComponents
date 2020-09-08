
import mechanical_components.meshes as meshes
import math

Rack1 = meshes.Rack(20/180.*math.pi)
input = {'z':13, 'db':40*1e-3, 'coefficient_profile_shift':0.3, 'rack':Rack1}
mesh1 = meshes.Mesh(**input)

d = mesh1.to_dict()
obj = meshes.Mesh.dict_to_object(d)

