use consts::*;
use model;
use raytracer;

use core::io::{Writer, WriterUtil};
use std;

fn write_ppm( fname: &str, width: uint, height: uint, pixels: &[raytracer::Color] ){

//    let writer = result::get( &io::file_writer( &Path(fname), [io::Create, io::Truncate] ) );
    let writer = result::get( &io::buffered_file_writer( &Path(fname) ) );
    writer.write_str(fmt!("P6\n%u %u\n255\n", width, height));
    for pixels.each |pixel| {
        writer.write([pixel.r, pixel.g, pixel.b]);
    };
}

#[main]
fn main()
{
    // Get command line args
    let args = os::args();

    if args.len() != 2u {
        io::println("Usage: rustray OBJ");
        io::println("");
        io::println("For example:");
        io::println("   $ wget http://groups.csail.mit.edu/graphics/classes/6.837/F03/models/cow-nonormals.obj");
        io::println("   $ ./rustray cow-nonormals.obj");
        io::println("   $ gimp oput.ppm");
        io::println("");
        fail!();
    }

    let start = std::time::precise_time_s();


    io::println(~"Reading \"" + args[1] + "\"...");
    let model = model::read_mesh( args[1] );
    
    let (depth,count) = model::count_kd_tree_nodes( &model.kd_tree );

    io::println(fmt!("Done.\nLoaded model.\n\tVerts: %u, Tris: %u\n\tKD-tree depth: %u, #nodes: %u",
                model.polys.vertices.len(),
                model.polys.indices.len()/3u,
                depth, count));

    io::print("Tracing rays... ");
    let start_tracing = std::time::precise_time_s();
    let pixels = raytracer::generate_raytraced_image(model, FOV, WIDTH, HEIGHT, SAMPLE_GRID_SIZE);
    io::println("Done!");
    let end_tracing = std::time::precise_time_s();
    
    let outputfile = "./oput.ppm";
    io::print(~"Writing \"" + outputfile + "\"...");
    write_ppm( outputfile, WIDTH, HEIGHT, pixels );
    io::println("Done!");

    let end = std::time::precise_time_s();
    io::print(fmt!("Total time: %3.3fs, of which tracing: %3.3f\n", end - start, end_tracing - start_tracing) );
}
