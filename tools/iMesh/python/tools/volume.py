from numpy import *
from numpy.linalg import *
from itaps import *
import sys
from pylab import *

def distance(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def tet_volume(coords):
    return abs(det( [coords[0]-coords[1],
                     coords[1]-coords[2],
                     coords[2]-coords[3]] )) / 6

def hex_volume(coords):
    # assumes not-quite logical vertex ordering
    def subvolume(a,b,c,d,e):
        base = ( distance(cross(b-a,d-a)) +
                 distance(cross(c-a,d-a)) ) / 2
        norm = cross(b-a,c-a)
        norm = norm / distance(norm)
        height = abs(dot(norm,e-a))
        return base*height / 3

    return subvolume(coords[0],coords[1],coords[3],coords[2],coords[7]) + \
           subvolume(coords[0],coords[1],coords[4],coords[5],coords[7]) + \
           subvolume(coords[1],coords[2],coords[5],coords[6],coords[7])

def calc_volume(filename):
    mesh = iMesh()
    mesh.load(mesh.rootSet, filename)

    volume = ndarray(mesh.getNumOfType(mesh.rootSet, iBase.type.region), float_)
    x=0
    for i in mesh.rootSet.iterate(iBase.type.region, iMesh.topology.all):
        topo = mesh.getEntTopo(i)
        curr = mesh.getVtxCoords( mesh.getEntAdj(i, iBase.type.vertex),
                                  iBase.storageOrder.interleaved )

        if topo == iMesh.topology.tetrahedron:
            volume[x] = tet_volume(curr)
        elif topo == iMesh.topology.hexahedron:
            volume[x] = hex_volume(curr)
        else:
            assert(False)
        x+=1
    return volume

if(len(sys.argv) != 3):
    print "Usage: python volume.py file1 file2"
    exit(1)

volume_pre  = calc_volume(sys.argv[1])
volume_post = calc_volume(sys.argv[2])
r = arange(len(volume_pre))

print volume_pre

plot(r, volume_pre,  linewidth=1)
plot(r, volume_post, linewidth=1)

xlabel('polyhedron index')
ylabel('volume')
title('Volume comparison pre- and post-deformation')
show()


