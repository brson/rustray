import io::writer_util;
import io::writer;
import raytracer::*;

const WIDTH : uint = 512u;
const HEIGHT : uint = 512u;
const FOV : float = 3.14159f / 3f ;

fn write_ppm( fname: str, width: uint, height: uint, pixels: [color] ){

	let writer = result::get( io::file_writer( fname, [io::create, io::truncate] ) );

	writer.write_str(#fmt("P6\n%u %u\n255\n", width, height));
	for pixel in pixels{
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
        io::println("   $ rustray cow-nonormals.obj");
        io::println("   $ gimp oput.ppm");
        io::println("");
        fail;
    }
	
    io::println("Reading \"" + args[1] + "\"...");
	let model = model::read_model( args[1] );
	let {depth,count} = model::count_kd_tree_nodes( model.kd_tree );

	io::println(#fmt("Done.\nLoaded model.\n\tVerts: %u, Tris: %u\n\tKD-tree depth: %u, #nodes: %u", 
				vec::len(model.mesh.vertices), 
				vec::len(model.mesh.indices)/3u,
				depth, count));

	io::print("Tracing rays... ");
	let pixels = raytracer::generate_raytraced_image(model, FOV, WIDTH, HEIGHT);
	io::println("Done!");
    
    let outputfile = "./oput.ppm";
	io::print("Writing \"" + outputfile + "\"...");
	write_ppm( outputfile, WIDTH, HEIGHT, pixels );	
	io::println("Done!");
}
