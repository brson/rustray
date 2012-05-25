rustray
January 2012

A raytracing proof-of-concept in Rust

DEMO:
   $ wget http://groups.csail.mit.edu/graphics/classes/6.837/F03/models/cow-nonormals.obj
   $ ./rustray cow-nonormals.obj
   Reading "cow-nonormals.obj"...
   Unrecognized line in .obj file: # The units used in this file are centimeters.
   Building kd-tree... Done.
   Loaded model.
       Verts: 4583, Tris: 5804
       KD-tree depth: 13, #nodes: 1905
   Tracing rays... Done!
   Writing "./oput.ppm"...Done!
   $ gimp oput.ppm


