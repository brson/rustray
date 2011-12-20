import std::{sort, io};
import math3d::*;

type mesh = { vertices: [vec3], indices: [uint], normals: [vec3] };

type model = { mesh: mesh, kd_tree: kd_tree };

fn find_split_plane( distances: [float], indices: [uint], faces: [uint] ) -> float {
	
	let face_distances = [];
	for f in faces {
		face_distances += [distances[indices[f*3u]]]; 
		face_distances += [distances[indices[f*3u+1u]]];
		face_distances += [distances[indices[f*3u+2u]]];
	}

	let sorted_distances = sort::merge_sort( {|a,b| a<b}, face_distances );
	sorted_distances[ vec::len(sorted_distances)/2u ]
}

fn split_triangles( splitter: float, distances: [float], indices: [uint], faces: [uint] ) -> {l:[uint],r:[uint]} {
	let l = [];
	let r = [];
	
	for f in faces {
		let d0 = distances[indices[f*3u   ]];
		let d1 = distances[indices[f*3u+1u]];
		let d2 = distances[indices[f*3u+2u]];

		let maxdist = float::max(d0, float::max(d1, d2));
		let mindist = float::min(d0, float::min(d1, d2));

		if mindist <= splitter {
			l += [f];
		}

		if maxdist > splitter {
			r += [f];
		}
	}

	{l:l, r:r}	
}

tag axis { 
	x;
	y;
	z; 
}


tag kd_tree {
	leaf( [uint] );
	node( axis, float, @kd_tree, @kd_tree );
}

fn build_kd_tree( 
	maxdepth: uint, 
	xdists: [float], 
	ydists: [float], 
	zdists: [float], 
	aabbmin: vec3,
	aabbmax: vec3,
	indices: [uint], 
	faces: [uint] ) -> @kd_tree {

	if maxdepth == 0u || vec::len(faces) <= 10u {
		ret @leaf( faces );
	}

	let extent = sub(aabbmax, aabbmin);
	let axis = extent.x > extent.y && extent.x > extent.z ? x :
		   extent.y > extent.z ? y : z;

	let dists;	
	alt axis {
		x() { dists = xdists }
		y() { dists = ydists }
		z() { dists = zdists }	
	};

	let s = find_split_plane( dists, indices, faces );	
	let {l,r} = split_triangles( s, dists, indices, faces );	

	// Stop when there's too much overlap between the two halves
	if vec::len(l) + vec::len(r) as float > vec::len(faces) as float * 1.35f {		
		ret @leaf( faces );
	}

	// adjust bounding boxes for children
	let left_aabbmax;
	let right_aabbmin;
	alt axis {
		x() { 
			left_aabbmax = { x: s with aabbmax };
			right_aabbmin = { x: s with aabbmin };
		}
		y() { 
			left_aabbmax = { y: s with aabbmax };
			right_aabbmin = { y: s with aabbmin };
		}
		z() { 
			left_aabbmax = { z: s with aabbmax };
			right_aabbmin = { z: s with aabbmin };
		}	
	};
	
	// build the two sub-trees
	let left_tree = build_kd_tree( 
				maxdepth - 1u, 
				xdists, 
				ydists, 
				zdists, 
				aabbmin, 
				left_aabbmax, 
				indices, 
				l );

	let right_tree = build_kd_tree( 
				maxdepth - 1u, 
				xdists, 
				ydists, 
				zdists, 
				right_aabbmin, 
				aabbmax, 
				indices, 
				r );


	@node( axis, s, left_tree, right_tree )
}

fn count_kd_tree_nodes( t: kd_tree ) -> {depth:uint, count:uint} {
	alt t {
		node(_,_,l,r) {
			let {depth:d0,count:c0} = count_kd_tree_nodes( *l );
			let {depth:d1,count:c1} = count_kd_tree_nodes( *r );
			ret {depth: float::max(d0, d1)+1u, count: c0+c1+1u };		
		}
		leaf {
			ret { depth:1u, count:1u };
		}
	}
}

fn read_model(fname: str) -> model {

	
	let m = read_mesh( fname );

	std::io::print("Building kd-tree... ");

	// just create a vector of 0..N-1 as our face array
	let max_tri_ix = vec::len(m.indices)/3u -1u;
	check( uint::le(0u, max_tri_ix) );
	let faces = vec::enum_uints(0u, max_tri_ix);

	// de-mux vertices for easier access later
	let xdists = [];
	let ydists = [];
	let zdists = [];
	let aabbmin = vec3(float::infinity, float::infinity, float::infinity);
	let aabbmax = vec3(float::neg_infinity, float::neg_infinity, float::neg_infinity);
	for v in m.vertices {
		xdists += [v.x];
		ydists += [v.y];
		zdists += [v.z];
		aabbmin = math3d::min(v, aabbmin);
		aabbmax = math3d::max(v, aabbmax);
	}

	let kdt = build_kd_tree( 
			20u, 
			xdists, 
			ydists, 
			zdists, 
			aabbmin,
			aabbmax,
			m.indices, 
			faces);
	
	{ mesh: m, kd_tree: *kdt }

}

fn read_mesh(fname: str) -> mesh {
	let reader = result::get( io::file_reader( fname ) );
	let vertices = [];
	let indices = [];
	
	let vert_normals : [mutable vec3]= [mutable];

	while !reader.eof() {
		let line : str = reader.read_line();
		if str::is_empty(line) {
			cont;
		}		
	
		let tokens = str::split(line, ' ' as u8 );	
		
		if tokens[0] == "v"{
			assert vec::len(tokens) == 4u;
			let v = vec3(	float::from_str(tokens[1]),
					float::from_str(tokens[2]),
					float::from_str(tokens[3]));
			assert v.x != float::NaN;
			assert v.y != float::NaN;
			assert v.z != float::NaN;

			vertices += [v];
			vert_normals += [mutable vec3(0f,0f,0f)];			

		} else if tokens[0] == "f" {
			assert vec::len(tokens) == 4u;
			let (i0,i1,i2) = (	uint::from_str(tokens[1])-1u,
						uint::from_str(tokens[2])-1u,
						uint::from_str(tokens[3])-1u );
			indices += [i0, i1, i2 ];
			let e1 = math3d::sub(vertices[i1], vertices[i0]);
			let e2 = math3d::sub(vertices[i2], vertices[i0]);
			let n = normalized(cross(e1,e2)); 
			
			vert_normals[i0] = math3d::add( vert_normals[i0], n );
			vert_normals[i1] = math3d::add( vert_normals[i1], n );
			vert_normals[i2] = math3d::add( vert_normals[i2], n );

			
		}
		else{
			io::println(#fmt("Unrecognized line in .obj file: %s", line));
		}
	}
	
	//for v in vert_normals {
 	//	v = normalized(v);
	//}

	ret {	vertices: vertices, 
		indices: indices, 
		normals: vec::map_mut( vert_normals, {|v| normalized(v)}) };
}



