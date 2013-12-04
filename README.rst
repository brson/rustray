Rustray
=======

May 2013

A raytracing proof-of-concept in Rust

Requirements
------------

- Rust-0.6 (incoming)

Compiling
---------
::

   $ rustc rustray.rc

It's usually worth compiling with '--opt-level 3'.

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
    Reading model file...Building kd-tree... Done.
    Loaded model.
	    Verts: 4583, Tris: 25811
	    KD-tree depth: 14, #nodes: 3341
    Tracing rays... using 4 tasks ... Done!
    Writing "./oput.ppm"...Done!
    Total time: 16.149s, of which tracing: 11.161
