#!/usr/bin/env python

# @file scale_stl.py
# 
# (c) Robert Bosch GmbH 2015. All rights reserved, also regarding any disposal, exploitation, reproduction, editing,  
# distribution, as well as in the event of applications for industrial property rights.
# 
# Scale mesh in stl file.
# 
# 
# @author Benjamin Maier (CR/AEG)
# 
# @date 05/2016
# 

import numpy as np
from numpy import sin, cos
import random
import sys
import csv
import struct

if __name__ == '__main__':
    
    # read in arguments
    if len(sys.argv) < 6:
      print "usage: {} <fileNameIn> <fileNameOut> [<x> <y> <z>]\n Scale vertices with factors [x,y,z].".format(sys.argv[0])
      quit()
    
    fileNameIn = sys.argv[1]
    fileNameOut = sys.argv[2]
    
    x = 0.0
    y = 0.0
    z = 0.0
    if len(sys.argv) == 6:
       x = float(sys.argv[3])
       y = float(sys.argv[4])
       z = float(sys.argv[5])
         
          
    print "\nConvert binary STL file {} to binary stl file {}\nScaling factors: [{},{},{}]".format(fileNameIn, fileNameOut, x, y, z)

    normalList = []
    vertexList = []

    # parse input file
    with open(fileNameIn, 'rb') as infile:
      
      header = infile.read(80)
      print "header: ",header
      
      nTriangles = struct.unpack('<i', infile.read(4))[0]
      print "nTriangles: ",nTriangles
            
      for i in range(nTriangles):
        normal = np.empty((3,1))
        normal[0] = struct.unpack('<f', infile.read(4))[0]
        normal[1] = struct.unpack('<f', infile.read(4))[0]
        normal[2] = struct.unpack('<f', infile.read(4))[0]
        vertex0 = np.empty((3,1))
        vertex1 = np.empty((3,1))
        vertex2 = np.empty((3,1))
        vertex0[0] = struct.unpack('<f', infile.read(4))[0]
        vertex0[1] = struct.unpack('<f', infile.read(4))[0]
        vertex0[2] = struct.unpack('<f', infile.read(4))[0]
        vertex1[0] = struct.unpack('<f', infile.read(4))[0]
        vertex1[1] = struct.unpack('<f', infile.read(4))[0]
        vertex1[2] = struct.unpack('<f', infile.read(4))[0]
        vertex2[0] = struct.unpack('<f', infile.read(4))[0]
        vertex2[1] = struct.unpack('<f', infile.read(4))[0]
        vertex2[2] = struct.unpack('<f', infile.read(4))[0]
        infile.read(2)
    
        normalList.append(normal)
        vertexList.append(vertex0)
        vertexList.append(vertex1)
        vertexList.append(vertex2)
      
      # compute center
      center = np.empty((3,1))
      for vertex in vertexList:
        
        center[0] = center[0]+vertex[0]
        center[1] = center[1]+vertex[1]
        center[2] = center[2]+vertex[2]
    
      print "number of vertices: ",len(vertexList)
    
      for i in range(3):
        center[i] /= len(vertexList)
      
      print "center of mass:",center
      
      # transform vertices
      for vertex in vertexList:
        vertex[0] *= x
        vertex[1] *= y
        vertex[2] *= z
        
      # compute bounding box
      minX = vertexList[0][0]
      maxX = vertexList[0][0]
      minY = vertexList[0][1]
      maxY = vertexList[0][1]
      minZ = vertexList[0][2]
      maxZ = vertexList[0][2]
      for vertex in vertexList:
        minX = min(minX, vertex[0])
        minY = min(minY, vertex[1])
        minZ = min(minZ, vertex[2])
        maxX = max(maxX, vertex[0])
        maxY = max(maxY, vertex[1])
        maxZ = max(maxZ, vertex[2])
      
      print "bounding box: [{},{}]x[{},{}]x[{},{}]".format(minX,maxX,minY,maxY,minZ,maxZ)
      
      # write to output file  
      outfile = open (fileNameOut, "wb")
      
      # header of binary file
      outfile.write("Binary STL file")

      for i in range(80-15):
        outfile.write(struct.pack("<b", 0))
      
      # number of triangles
      outfile.write(struct.pack("<I", nTriangles))
            
      # loop over points
      for i in range(nTriangles):
        
        a = np.array([vertexList[3*i+0][0][0], vertexList[3*i+0][1][0], vertexList[3*i+0][2][0]])
        b = np.array([vertexList[3*i+1][0][0], vertexList[3*i+1][1][0], vertexList[3*i+1][2][0]])
        c = np.array([vertexList[3*i+2][0][0], vertexList[3*i+2][1][0], vertexList[3*i+2][2][0]])
        
        # compute normal
        n = np.cross(-a+b, -a+c)
        n = n / np.linalg.norm(n)
        
        outfile.write(struct.pack("<f", n[0]))
        outfile.write(struct.pack("<f", n[1]))
        outfile.write(struct.pack("<f", n[2]))
                
        outfile.write(struct.pack("<f", a[0]))
        outfile.write(struct.pack("<f", a[1]))
        outfile.write(struct.pack("<f", a[2]))
        
        outfile.write(struct.pack("<f", b[0]))
        outfile.write(struct.pack("<f", b[1]))
        outfile.write(struct.pack("<f", b[2]))
        
        outfile.write(struct.pack("<f", c[0]))
        outfile.write(struct.pack("<f", c[1]))
        outfile.write(struct.pack("<f", c[2]))
        
        outfile.write(struct.pack("<h", 0))
        
      # close binary file
      outfile.close()
      
      print "File {} written.".format(fileNameOut)
        
      print "number of triangles: {}".format(nTriangles)
