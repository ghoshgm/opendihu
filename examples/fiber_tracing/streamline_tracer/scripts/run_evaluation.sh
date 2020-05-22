#!/bin/bash
#
# This scripts tests the different algorithms for create_mesh.py

# directories
basedir=$(pwd)/..
opendihu_directory=$(pwd)/../../../..
stl_utility_directory=${opendihu_directory}/scripts/stl_utility
pyod=${opendihu_directory}/dependencies/python/install/bin/python3

# ------------
# settings

# z range of the mesh for which fibers should be extracted
min_z=40
max_z=260

# number of subdivisions of the z range, this has to be odd for quadratic elements
#n_rings=13
#n_rings=101
n_rings=43

# number of fibers in x or y direction,  has to be even for quadratic elements
#n_points_x=100
#n_points_x=32
n_points_x=10
#n_points_x=4

# input mesh to use, this is an STL mesh of the surface of the geometry
input_file=${basedir}/original_meshes/biceps_full.stl
# -----------------

# change to working directory for meshes
mkdir -p ${basedir}/processed_meshes
cd ${basedir}/processed_meshes

cp $input_file mesh_00.stl    # cp input file

# ---- python part, preprocessing ------

# cut surface mesh at specified z positions and create rings from it. Write result to `rings_created`, debugging output to out/mesh_01.stl (takes <1min)
# usage: ./create_rings.py <input file> [<n rings> [<min z> <max z>]]
$pyod ${basedir}/scripts/utility/create_rings.py $input_file $n_rings $min_z $max_z

# extract the existing rings from the surface mesh, i.e. use the nodes in the STL mesh
# (this is the alternative to create_rings.py)   
# write result to `rings_extracted`.
#$pyod ${basedir}/scripts/utility/extract_rings.py ../biceps.stl     

# rename ring output file to `rings`
mv rings_created rings      
#mv rings_extracted rings

#######################################################################################################
# loop over triangulation type
for triangulation_type in 2; do  # triangulation_type:  0 = scipy, 1 = triangle, 2 = custom (1 is best)

  # loop over parametric space shape
  for parametric_space_shape in 4; do   # 0 = unit circle, 1 = unit square, 2 = unit square with adjusted grid, 3 = unit circle with adjusted grid

    # create a mesh, reads input from `rings`, write output to `mesh`, debugging output to out/mesh_02* - out/mesh_09* (takes ~2min)
    # usage: ./create_mesh.py [<triangulation_type> [<parametric_space_shape> [<n_points_x> [<n_grid_points_x>]]]]"

    $pyod ${basedir}/scripts/utility/create_mesh.py $triangulation_type $parametric_space_shape $n_points_x

    cp -r out out_$triangulation_type_$parametric_space_shape

  done

done
