import re
import os
import sys
import subprocess
import csv
from datetime import *
sys.path.append('../build/lib.linux-x86_64-2.5/')

from itaps import *

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


c_stats  = {}
py_stats = {}
list_stats = {}

path = raw_input('Location of test files: ')
if(len(path) == 0):
    path = '.'

##### Run some tests in C #####
lines = csv.reader( subprocess.Popen('./perf %s' % path, shell=True,
                                     stdout=subprocess.PIPE).stdout )
for row in lines:
    c_stats[row[0]] = [float(i)/100 for i in row[1:]]


##### Run some tests in Python #####
os.chdir(path)
timer = stopwatch()

for file in filter(testable, os.listdir(path)):
    py_stats[file] = []
    list_stats[file] = []

    ##### 1 #####
    timer.reset()
    for x in range(100):
        m = iMesh()
        m.load(m.rootSet, file)
        m = None
    py_stats[file].append( timer.delta()/100 )
    list_stats[file].append(0)

    ##### Intermission #####
    mesh = iMesh()
    root = mesh.rootSet
    mesh.load(root, file)

    ##### 2 #####
    timer.reset()
    for x in range(100):
        mesh.getAdjEntIndices(root, iBase.type.all,
                              iMesh.topology.all, iBase.type.all)
    py_stats[file].append( timer.delta()/100 )
    list_stats[file].append(0)

    ##### Intermission #####
    arr = mesh.getAdjEntIndices(root, iBase.type.all,
                                iMesh.topology.all, iBase.type.all)
    list = map(lambda x: x.tolist(), arr)

    ##### 3 #####
    timer.reset()
    for x in range(100):
        mesh.getEntAdj(arr[0], iBase.type.all)
    py_stats[file].append( timer.delta()/100 )

    timer.reset()
    for x in range(100):
        mesh.getEntAdj(list[0], iBase.type.all)
    list_stats[file].append( timer.delta()/100 )

    ##### 4 #####
    timer.reset()
    for x in range(100):
        mesh.getEnt2ndAdj(arr[0], iBase.type.edge, iBase.type.vertex)
    py_stats[file].append( timer.delta()/100 )

    timer.reset()
    for x in range(100):
        mesh.getEnt2ndAdj(list[0], iBase.type.edge, iBase.type.vertex)
    list_stats[file].append( timer.delta()/100 )

    mesh = None
    arr = None
    list = None


# for some reason, importing this earlier makes pytaps run faster????
import matplotlib.pyplot as plt
import numpy

i = 1
for file in filter(testable, os.listdir(path)):
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
