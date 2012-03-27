import std::*;
import float::{fmax,fmin};
import math3d::*;

type color = { r:u8, g:u8, b:u8 };

fn for_each_pixel( width: uint, height: uint, f : fn (x: uint, y: uint) -> color ) -> [color]{
	let mut img_pixels = [];


	uint::range( 0u, height ) { |row|	
		uint::range(0u, width) { |column|
			img_pixels += [f(column,row)];				
		}
	}
	ret img_pixels;
}

fn generate_test_image( width: uint, height: uint ) -> [color]{
	for_each_pixel( width, height, {|i,j|
		let checker = (i / 32u + j / 32u) % 2u == 0u;
		let ifloat = (i as float) * 255.0f / (width as float);
		let jfloat = (j as float) * 255.0f / (height as float);
		if checker { 
			{r: 0u8, g: 0u8, b:155u8} 
		} else { 
			{r: ifloat as u8, 
			 g: jfloat as u8,
			 b:0u8} 
		}
	})
}
 

fn get_ray( horizontalFOV: float, width: uint, height: uint, x: uint, y: uint) -> ray{
	let dirx = (x as float) - ((width/2u) as float);
	let diry = -((y as float) - ((height/2u) as float));
	let dirz = -((width/2u) as float) / float::tan(horizontalFOV*0.5);
	{ origin: vec3(0.5f, 0f, 10f), 
	  dir: normalized( vec3( dirx+0.00001f, 
				 diry+0.00001f, 
				 dirz) ) }
}

fn get_triangle( m : model::mesh, ix : uint ) -> triangle{
	{ 	p1: m.vertices[ m.indices[ix*3u   ] ],
		p2: m.vertices[ m.indices[ix*3u+1u] ],
		p3: m.vertices[ m.indices[ix*3u+2u] ] }
}

fn clamp( x: float, lo: float, hi: float ) -> float{
	fmin(hi, fmax( x, lo ))
}

fn abs( f:float )->float{
        if f < 0f { -f } else { f }
}

fn trace_kd_tree( 
	mesh: model::mesh, 
	kd_tree: model::kd_tree, 
	r: ray, 
	mint: float, 
	maxt: float ) 
	-> option<(hit_result, uint)> {

	alt kd_tree {
		model::leaf(tris) {

			let mut res : option<(hit_result, uint)> = option::none;
			for tri_index in tris{
				let t = get_triangle( mesh, tri_index ); 
				let new_hit_result = ray_triangle_intersect( r, t );

				alt (res, new_hit_result){
					(option::none(), option::some(hr)) {
						res = option::some((hr,tri_index));
					}
					(option::some((hr1,_)), option::some(hr2)) {
						if hr1.t > hr2.t {
							res = option::some((hr2,tri_index));
						}
					}
					_ {}
				}
			}

			ret res;
		}
		model::node(axis, splitter, left_tree, right_tree) {
				
			// find the scalar direction/origin for the current axis
			let (dir, origin) = alt axis {
				model::x() { (r.dir.x, r.origin.x) }
				model::y() { (r.dir.y, r.origin.y) }
				model::z() { (r.dir.z, r.origin.z) }
			};

			// find intersection with plane
			// origin + dir*t = splitter

			let plane_dist = (splitter - origin) / dir;
							
	        let (near,far) = if origin < splitter || (origin == plane_dist && dir <= 0f) {
 		        (left_tree,right_tree)
            } else {
                (right_tree, left_tree)
            };

			if plane_dist > maxt || plane_dist < 0f {
				ret trace_kd_tree( mesh, *near, r, mint, maxt);
			} else if plane_dist < mint {
				ret trace_kd_tree(mesh, *far, r, mint, maxt);
			} else{
				let near_hit = trace_kd_tree( mesh, *near, r, mint, plane_dist);
				
				alt near_hit {
					option::some((h1,_)) {

						// if we hit something before the splitting plane, 
						// we can early out now.
						if h1.t < plane_dist { 
							ret near_hit; 
						}
						
						// else, check far node too...
						let far_hit = trace_kd_tree(mesh, *far, r, plane_dist, maxt);
						alt far_hit {
						  option::some((h2,_)) if h2.t < h1.t {
								ret far_hit;
							}
						  _ { ret near_hit; }
						}							
					}
					_ { 
						ret trace_kd_tree(mesh, *far, r, plane_dist, maxt); 
					}
				}	
			}
					
							
		}
	}
}

fn trace_soup( mesh: model::mesh, r: ray) -> option<(hit_result, uint)>{
	
	let mut res : option<(hit_result, uint)> = option::none;
	
	uint::range( 0u, vec::len( mesh.indices ) / 3u) { |tri_ix|
		let tri = get_triangle(mesh,tri_ix);

		let new_hit = ray_triangle_intersect( r, tri );

		alt (res, new_hit) {
			(option::none(),option::some(hit)) { 
				res = option::some((hit, tri_ix)); 
			}
			(option::some((old_hit,_)), option::some(hit))
		                if hit.t < old_hit.t && hit.t > 0f {
				res = option::some((hit, tri_ix));
			}
			_ {}
		}						
	} 

	ret res;
}

fn shade( pos: vec3, n: vec3, r: ray, color: vec3, reflectivity: float, 
	  shadow_test: fn( light_dir: vec3 ) -> bool, 
	  reflection_color: fn( reflection_dir: vec3) -> vec3 ) -> vec3 {

	let light_pos : vec3 = {x: 3f, y: 3f, z: 8.5f};
	let light_strength = 50f;

	let light_vec = sub(light_pos, pos);
	let view_vec = normalized(sub(r.origin, pos));

	let direct_light = if shadow_test(light_vec){
		vec3(0f, 0f, 0f)
	} else {

		let light_vec_n = normalized(light_vec);

		
		let half_vector = normalized(add(light_vec_n, view_vec));

		let s = dot(n,half_vector);
		let specular = s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s*s;				

	
		let atten = light_strength*(1.0f/length_sq(light_vec) + specular*0.01f);
	
		let intensity = atten*dot(n, light_vec_n );
	
		scale(color, intensity)
	};

	let reflection = add(scale(n, dot(view_vec, n)*2f), view_vec);
	let rcolor = if reflectivity > 0.01f { reflection_color( reflection ) } else { vec3(0f,0f,0f) };

	lerp(direct_light, rcolor, reflectivity)
}


type intersection = { pos: vec3, n: vec3, color: vec3, reflectivity: float };

fn trace_ray( r : ray, model : model::model, mint: float, maxt: float) -> option<intersection> {

	
	let use_kd_tree = true;	

	// trace against checkerboard first
	let checkerboard_height = -3.5f;
	let checker_hit_t = (checkerboard_height - r.origin.y) / r.dir.y;

	// compute checkerboard color, if we hit the floor plane
	let mut (result, new_maxt) = if checker_hit_t > mint && checker_hit_t < maxt {

			let s = 1.0f;			
			let pos = add(r.origin, scale(r.dir, checker_hit_t));
			
			// hacky checkerboard pattern
			let (u,v) = ((pos.x*s + 100000f) as uint, (pos.z*s  + 100000f) as uint);
        	let color = if (u + v) % 2u == 0u { vec3(1f,1f,1f) } else { vec3(0.5f,0.5f,0.2f) };
			let intersection = option::some( { 
						pos: pos,
						n: vec3(0f,1f,0f),
						color: color,
						reflectivity: 0.0f } );
			(intersection, checker_hit_t)
		} else {
			(option::none, maxt)
		};
	
	// trace against scene
	let trace_result = if use_kd_tree {
				trace_kd_tree( model.mesh, model.kd_tree, r, mint, new_maxt )
    } else {
				trace_soup( model.mesh, r)
    };

	option::may( trace_result, { |hit|
		let (hit_info, tri_ix) = hit;
	
	
		if hit_info.t > 0f {
			let pos = add( r.origin, scale(r.dir, hit_info.t));
						
			let (v0,v1,v2) = (	model.mesh.indices[tri_ix*3u], 
						model.mesh.indices[tri_ix*3u+1u], 	
						model.mesh.indices[tri_ix*3u+2u] );		

			// interpolate vertex normals...
			let n = normalized(
					add( scale( model.mesh.normals[v0], hit_info.barycentric.z),
					add( scale( model.mesh.normals[v1], hit_info.barycentric.x),
				     	     scale( model.mesh.normals[v2], hit_info.barycentric.y))));			
					

			result = option::some( {
					pos: pos,
					n: n,
					color: vec3(1f,1f,1f),
					reflectivity: 0.3f } );	

		}

	}); 


	ret result;
}

fn get_color( r: ray, model: model::model, tmin: float, tmax: float, max_depth: uint) -> vec3 {
	let default_color = vec3(0.25f, 0.25f, 0.45f);	
	
	if max_depth == 0u {
		ret default_color;
	}
	
	alt trace_ray( r, model, tmin, tmax ) {
		option::some({pos,n,color, reflectivity}) { 
			let surface_origin = add(pos, scale(n,0.00001f));

			shade(pos, n, r, color, reflectivity, {|light_vec| 
				let shadow_ray = {origin: surface_origin, dir: light_vec};
				let st = trace_ray(shadow_ray, model, 0f, 1f);
				option::is_some(st)
			}, {|reflection_dir| 
				let reflection_ray = {origin: surface_origin, dir: reflection_dir};
				get_color(reflection_ray, model, tmin, tmax, max_depth - 1u)
			})
			
		}
		_ { default_color }
	}

}

fn generate_raytraced_image( 
	model: model::model, 
	horizontalFOV: float, 
	width: uint, 
	height: uint ) -> [color] 
{	
	for_each_pixel( width, height, {|x,y|
		let r = get_ray(horizontalFOV, width, height, x,y);

		let shaded_color = scale(get_color(r, model, 0f, float::infinity, 5u), 255f);

		{	r: clamp(shaded_color.x, 0f, 255f) as u8, 
			g: clamp(shaded_color.y, 0f, 255f) as u8, 
			b: clamp(shaded_color.z, 0f, 255f) as u8 }
	})
}

