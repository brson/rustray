Rustray
=======

October 2012

A raytracing proof-of-concept in Rust

Requirements
------------

- Rust-0.4

Compiling
---------
::

   $ rustc rustray.rc

Quality Settings
----------------

Now raytracing parameters are hard-coded in ``consts.rs``,
just edit it and re-compile ``rustray`` again.

Demo
----
This might take some time even on the strong machine, so consider choosing more
simple model or adjust quality settings.
::

   $ wget http://groups.csail.mit.edu/graphics/classes/6.837/F03/models/cow-nonormals.obj
   $ ./rustray cow-nonormals.obj
   Reading "cow-nonormals.obj"...
   Building kd-tree... Done.
   Loaded model.
       Verts: 4583, Tris: 5804
       KD-tree depth: 13, #nodes: 1905
   Tracing rays... Done!
   Writing "./oput.ppm"...Done!
   $ gimp oput.ppm
