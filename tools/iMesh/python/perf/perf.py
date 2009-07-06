import re
import os
import subprocess
from datetime import *
from itaps import *
from optparse import OptionParser

class stopwatch:
    def __init__(self):
        self.reset()

    def delta(self):
        d = datetime.now() - self.start
        return d.seconds + d.microseconds/1000000.0

    def reset(self):
        self.start = datetime.now()

def testable(name):
    return re.search(r'\.g$',name) != None and \
        not re.match('testquad',name)


##### Parse command line options #####
parser = OptionParser(usage='usage: %prog [options] [files...]')
parser.add_option('-c', '--count', action='store', type='int', default=20,
                  help='number of times to run the performance tests',
                  metavar='NUM')

(options, args) = parser.parse_args()

if len(args) == 0:
    args = ['.']

files = []
for path in args:
    if os.path.isdir(path):
        files.extend([ os.path.join(path,i) for i in 
                       filter(testable, os.listdir(path)) ])
    else:
        files.append(path)

count = options.count


c_stats    = {}
py_stats   = {}
list_stats = {}

##### Run some tests in C #####
for file in files:
    results = subprocess.Popen('./perf %d "%s"' % (count, file), shell=True,
                               stdout=subprocess.PIPE).stdout
    c_stats[file] = [float(line)/count for line in results]

##### Run some tests in Python #####
timer = stopwatch()

for file in files:
    py_stats[file] = []
    list_stats[file] = []

    ##### 1 #####
    timer.reset()
    for x in range(count):
        m = iMesh.Mesh()
        m.load(file)
        m = None
    py_stats[file].append( timer.delta()/count )
    list_stats[file].append(0)

    ##### Intermission #####
    mesh = iMesh.Mesh()
    mesh.load(file)

    ##### 2 #####
    timer.reset()
    for x in range(count):
        mesh.getAdjEntIndices(iBase.Type.all, iMesh.Topology.all,
                              iBase.Type.all)
    py_stats[file].append( timer.delta()/count )
    list_stats[file].append(0)

    ##### Intermission #####
    arr = mesh.getEntities(iBase.Type.all, iMesh.Topology.all)
    list = arr.tolist()

    ##### 3 #####
    timer.reset()
    for x in range(count):
        mesh.getEntAdj(arr, iBase.Type.all)
    py_stats[file].append( timer.delta()/count )

    timer.reset()
    for x in range(count):
        mesh.getEntAdj(list, iBase.Type.all)
    list_stats[file].append( timer.delta()/count )

    ##### 4 #####
    timer.reset()
    for x in range(count):
        mesh.getEnt2ndAdj(arr, iBase.Type.edge, iBase.Type.vertex)
    py_stats[file].append( timer.delta()/count )

    timer.reset()
    for x in range(count):
        mesh.getEnt2ndAdj(list, iBase.Type.edge, iBase.Type.vertex)
    list_stats[file].append( timer.delta()/count )

    mesh = None
    arr  = None
    list = None


# for some reason, importing this earlier makes pytaps run faster????
import matplotlib.pyplot as plt
import numpy

i = 1
for file in files:
    ind = numpy.arange(len(py_stats[file]))
    width = 0.25
    plt.figure(i)
    plt.title(file)
    plt.ylabel('time (s)')
    plt.xticks(ind+1.5*width,
               ('load()','getAdjEntIndices()','getEntAdj()','getEnt2ndAdj()',
                'getEntities()'))

    bars1 = plt.bar(ind,         py_stats[file],   width, color='r')
    bars2 = plt.bar(ind+width,   c_stats[file],    width, color='b')
    bars3 = plt.bar(ind+2*width, list_stats[file], width, color='g')

    plt.legend( (bars1[0], bars2[0], bars3[0]),
                ('Python (array)', 'C', 'Python (list)'), loc=0 )

    i+=1

plt.show()
