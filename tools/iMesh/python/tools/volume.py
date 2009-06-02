from numpy import *
from numpy.linalg import *
from itaps import *
from optparse import OptionParser
from pylab import *

def distance(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def tet_volume(coords):
    return abs(det( [coords[0]-coords[1],
                     coords[1]-coords[2],
                     coords[2]-coords[3]] )) / 6

def hex_volume(coords):
    # assumes not-quite logical vertex ordering
    def subvolume(a, b, c, d, e):
        base = ( distance(cross(b-a, d-a)) +
                 distance(cross(c-a, d-a)) ) / 2
        norm = cross(b-a, c-a)
        norm = norm / distance(norm)
        height = abs(dot(norm, e-a))
        return base*height / 3

    return subvolume(coords[0], coords[1], coords[3], coords[2], coords[7]) + \
           subvolume(coords[0], coords[1], coords[4], coords[5], coords[7]) + \
           subvolume(coords[1], coords[2], coords[5], coords[6], coords[7])

def calc_volume(mesh):
    volume = ndarray(mesh.rootSet.getNumOfType(iBase.Type.region), float_)
    x=0
    for i in mesh.rootSet.iterate(iBase.Type.region, iMesh.Topology.all):
        topo = mesh.getEntTopo(i)
        curr = mesh.getVtxCoords( mesh.getEntAdj(i, iBase.Type.vertex),
                                  iBase.StorageOrder.interleaved )

        if topo == iMesh.Topology.tetrahedron:
            volume[x] = tet_volume(curr)
        elif topo == iMesh.Topology.hexahedron:
            volume[x] = hex_volume(curr)
        else:
            assert(False)
        x+=1
    return volume

parser = OptionParser()
parser.add_option('-r', '--raw', action='store_true', dest='raw')

(options, args) = parser.parse_args()

if len(args) != 2:
    print 'Usage: python volume.py [opts] file1 file2'
    exit(1)

mesh_pre = iMesh.Mesh()
mesh_pre.load(args[0])
mesh_post = iMesh.Mesh()
mesh_post.load(args[1])

if mesh_pre. rootSet.getNumOfType(iBase.Type.region) != \
   mesh_post.rootSet.getNumOfType(iBase.Type.region):
    print 'volume.py: Meshes should have the same number of regions'
    exit(1)

volume_pre  = calc_volume(mesh_pre)
volume_post = calc_volume(mesh_post)

volume_diff = volume_pre / volume_post

if options.raw:
    for i in range(len(volume_pre)):
        print '%f,%f,%f' % (volume_pre[i], volume_post[i], volume_diff[i])
else:
    r = arange(len(volume_pre))

    subplot(2,1,1)
    plot(r, volume_pre,  '.',
         r, volume_post, '.')

    title('Volume comparison pre- and post-deformation')

    xlabel('polyhedron index')
    ylabel('volume')

    subplot(2,1,2)
    plot(r, volume_diff, '.')

    xlabel('polyhedron index')
    ylabel('volume ratio')

    show()
