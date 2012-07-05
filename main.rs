import io::writer_util;
import io::writer;
import raytracer::*;

// TODO: These things should all be overridable by the command line
const WIDTH : uint = 256u;
const HEIGHT : uint = 256u;
const FOV : f32 = 3.14159f32 / 3f32 ;
const SAMPLE_GRID_SIZE : uint = 1u;
const NUM_GI_SAMPLES_SQRT: uint = 8u;
const NUM_LIGHT_SAMPLES : uint = 32u;
const MAX_TRACE_DEPTH : uint = 3u;
const USE_SMOOTH_NORMALS_FOR_GI : bool = true;
const USE_SMOOTH_NORMALS_FOR_DIRECT_LIGHTING : bool = true;

fn write_ppm( fname: str, width: uint, height: uint, pixels: [color] ){

	let writer = result::get( io::file_writer( fname, [io::create, io::truncate] ) );

	writer.write_str(#fmt("P6\n%u %u\n255\n", width, height));
	for pixels.each |pixel| {
		writer.write([pixel.r, pixel.g, pixel.b]);
	};
}

fn main( args: [str] )
{
    if vec::len(args) != 2u {
        io::println("Usage: rustray OBJ");
        io::println("");
        io::println("For example:");
        io::println("   $ wget http://groups.csail.mit.edu/graphics/classes/6.837/F03/models/cow-nonormals.obj");
        io::println("   $ ./rustray cow-nonormals.obj");
        io::println("   $ gimp oput.ppm");
        io::println("");
        fail;
    }
	
	let start = std::time::precise_time_s();

	
    io::println("Reading \"" + args[1] + "\"...");
	let model = model::read_mesh( args[1] );
	let {depth,count} = model::count_kd_tree_nodes( model.kd_tree );

	io::println(#fmt("Done.\nLoaded model.\n\tVerts: %u, Tris: %u\n\tKD-tree depth: %u, #nodes: %u", 
				vec::len(model.polys.vertices), 
				vec::len(model.polys.indices)/3u,
				depth, count));

	io::print("Tracing rays... ");
	let start_tracing = std::time::precise_time_s();
	let pixels = raytracer::generate_raytraced_image(model, FOV, WIDTH, HEIGHT, SAMPLE_GRID_SIZE);
	io::println("Done!");
    
    let outputfile = "./oput.ppm";
	io::print("Writing \"" + outputfile + "\"...");
	write_ppm( outputfile, WIDTH, HEIGHT, pixels );	
	io::println("Done!");
	
	let end = std::time::precise_time_s();
    io::print(#fmt("Total time: %3.3fs, of which tracing: %3.3f\n", end - start, end - start_tracing) );
}
