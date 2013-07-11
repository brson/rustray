use consts::*; // for the consts.. ugh... make the constants go away
use math3d::*;
use concurrent;
use model;
use std::*;
use std::rand::{RngUtil,Rng,task_rng};

use extra;

pub struct Color { r:u8, g:u8, b:u8 }

#[inline]
fn for_each_pixel( width: uint, height: uint, f : &fn (x: uint, y: uint) -> Color ) -> ~[Color]{
    let mut img_pixels = vec::with_capacity(height*width);

    for uint::range( 0, height ) |row| {
        for uint::range(0, width) |column| {
            img_pixels.push(f(column,row));
        }
    }
    img_pixels
}

#[inline]
fn get_ray( horizontalFOV: f32, width: uint, height: uint, x: uint, y: uint, sample_jitter : (f32,f32)) -> Ray {
    let (jitterx,jittery) = sample_jitter;
    let dirx = (x as f32) - ((width/2) as f32) + jitterx;
    let diry = -((y as f32) - ((height/2) as f32)) + jittery;
    let dirz = -((width/2u) as f32) / (horizontalFOV*0.5).tan();
    Ray{ origin: vec3(0.0, 0.0, 1.0),
      dir: normalized( vec3( dirx, diry, dirz) ) }
}

#[deriving(Clone)]
struct rand_env {
    floats: ~[f32],
    disk_samples: ~[(f32,f32)],
    hemicos_samples: ~[vec3]
}

fn get_rand_env() -> rand_env {
    let mut gen = task_rng();

    let disk_samples = do vec::from_fn(513u) |_x| {
        // compute random position on light disk
        let r_sqrt = gen.gen::<f32>().sqrt();
        let theta = gen.gen::<f32>() * 2.0 * f32::consts::pi;
        (r_sqrt * theta.cos(), r_sqrt*theta.sin())
    };

    let mut hemicos_samples = vec::with_capacity(NUM_GI_SAMPLES_SQRT * NUM_GI_SAMPLES_SQRT);

    for uint::range(0, NUM_GI_SAMPLES_SQRT) |x| {
        for uint::range(0, NUM_GI_SAMPLES_SQRT) |y| {
            let (u,v) = (    ( (x as f32) + gen.gen() ) / (NUM_GI_SAMPLES_SQRT as f32),
                            ( (y as f32) + gen.gen() ) / (NUM_GI_SAMPLES_SQRT as f32) );
            hemicos_samples.push(cosine_hemisphere_sample(u,v));
        }
    };

    rand_env{
        floats: vec::from_fn(513u, |_x| gen.gen() ),
        disk_samples: disk_samples,
        hemicos_samples: hemicos_samples }
}

#[inline]
fn sample_floats_2d_offset( offset: uint, rnd: &rand_env, num: uint, body: &fn(f32,f32) ) {
    let mut ix = offset % rnd.floats.len();
    for num.times {
        let r1 = rnd.floats[ix];
        ix = (ix + 1) % rnd.floats.len();
        let r2 = rnd.floats[ix];
        body(r1,r2);
        ix = (ix + 1) % rnd.floats.len();
    }
}

#[inline]
fn sample_disk( rnd: &rand_env, num: uint, body: &fn(f32,f32) ){
    let mut rng = rand::task_rng();
    if ( num == 1 ) {
        body(0.0,0.0);
    } else {
        let mut ix = rng.gen::<uint>() % rnd.disk_samples.len(); // start at random location
        for num.times {
            let (u,v) = rnd.disk_samples[ix];
            body(u,v);
            ix = (ix + 1) % rnd.disk_samples.len();
        };
    }
}

// Sample 2 floats at a time, starting with an offset that's passed in
#[inline]
fn sample_floats_2d( rnd: &rand_env, num: uint, body: &fn(f32,f32) ) {
    let mut rng = rand::task_rng();
    sample_floats_2d_offset( rng.next() as uint, rnd, num, body);
}

#[inline(always)]
fn sample_stratified_2d( rnd: &rand_env, m: uint, n : uint, body: &fn(f32,f32) ) {
    let m_inv = 1.0/(m as f32);
    let n_inv = 1.0/(n as f32);
    let mut rng = rand::task_rng();
    let start_offset = rng.next();
    for uint::range( 0, m ) |samplex| {
        // sample one "row" of 2d floats
        let mut sampley = 0;
        do sample_floats_2d_offset( (start_offset as uint) + (n*samplex as uint), rnd, n ) |u,v| {
            body(  ((samplex as f32) + u) * m_inv,
                    ((sampley as f32) + v) * n_inv );
            sampley += 1;
        };
    }
}

#[inline(always)]
fn sample_cosine_hemisphere( rnd: &rand_env, n: vec3, body: &fn(vec3) ) {
    let mut rng = rand::task_rng();
    let rot_to_up = rotate_to_up(n);
    let random_rot = rotate_y( rnd.floats[ rng.next() as uint % rnd.floats.len() ] ); // random angle about y
    let m = mul_mtx33(rot_to_up, random_rot);
        for rnd.hemicos_samples.iter().advance |s| {
        body(transform(m,*s));
    }
}

#[inline(always)]
fn get_triangle( m : &model::polysoup, ix : uint ) -> Triangle{
    Triangle{   p1: m.vertices[ m.indices[ix*3   ] ],
                p2: m.vertices[ m.indices[ix*3+1] ],
                p3: m.vertices[ m.indices[ix*3+2] ] }
}

#[inline(always)]
fn clamp( x: f32, lo: f32, hi: f32 ) -> f32{
    //f32::min(hi, f32::max( x, lo ))
    if x < lo { lo } else if x > hi { hi } else { x }
}

#[inline]
fn trace_kd_tree(
    polys: &model::polysoup,
    kd_tree_nodes: &[model::kd_tree_node],
    kd_tree_root: uint,
    r: &Ray,
    inv_dir: vec3,
    inmint: f32,
    inmaxt: f32 )
    -> Option<(HitResult, uint)> {

    let mut res : Option<(HitResult, uint)> = option::None;
    let mut closest_hit = inmaxt;

    let mut stack : ~[(uint, f32, f32)] = ~[];
    let mut mint = inmint;
    let mut maxt = inmaxt;
    let mut cur_node = kd_tree_root;

    loop {
        // skip any nodes that have been superceded
        // by a closer hit.
        while mint >= closest_hit {
            if ( stack.len() > 0 ){
                let (n,mn,mx) = stack.pop();
                cur_node = n;
                mint = mn;
                maxt = mx;
            } else {
                return res;
            }
        }

        match kd_tree_nodes[cur_node] {
            model::leaf(tri_begin, tri_count) => {
                let mut tri_index : u32 = tri_begin;
                while tri_index < tri_begin+tri_count {

                    let t = &get_triangle( polys, tri_index as uint );
                    let new_hit_result = r.intersect(t);

                    match (res, new_hit_result){
                        (option::None(), option::Some(hr)) => {
                            res = option::Some((hr,tri_index as uint));
                            closest_hit = hr.t;
                        }
                        (option::Some((hr1,_)), option::Some(hr2)) if hr1.t > hr2.t => {
                            res = option::Some((hr2,tri_index as uint));
                            closest_hit = hr2.t;
                        }
                        _ => {}
                    }
                    tri_index += 1;
                }

                if ( stack.len() > 0 ){
                    let (n,mn,mx) = stack.pop();
                    cur_node = n;
                    mint = mn;
                    maxt = mx;
                } else {
                    return res;
                }
            }
            model::node(axis, splitter, right_tree) => {
                // find the scalar direction/origin for the current axis
                let (inv_dir_scalar, origin) = match axis {
                    model::x() => { (inv_dir.x, r.origin.x) }
                    model::y() => { (inv_dir.y, r.origin.y) }
                    model::z() => { (inv_dir.z, r.origin.z) }
                };
                // figure out which side of the spliting plane the ray origin is
                // i.e. which child we need to test first.
                let (near,far) = if origin < splitter || (origin == splitter && inv_dir_scalar >= 0.0) {
                    ((cur_node+1) as uint,right_tree as uint)
                } else {
                    (right_tree as uint, (cur_node+1) as uint)
                };
                // find intersection with plane
                // origin + dir*plane_dist = splitter
                let plane_dist = (splitter - origin) * inv_dir_scalar;

                if plane_dist > maxt || plane_dist <= 0.0 {
                    cur_node = near;
                } else if plane_dist < mint {
                    cur_node = far;
                } else{
                    stack.push((far, plane_dist, maxt) );
                    cur_node = near;
                    maxt = plane_dist;
                }
            }
        }
    }
}

#[inline]
fn trace_kd_tree_shadow(
    polys: &model::polysoup,
    kd_tree_nodes: &[model::kd_tree_node],
    kd_tree_root: uint,
    r: &Ray,
    inv_dir: vec3,
    inmint: f32,
    inmaxt: f32 )
    -> bool {

    let mut stack : ~[(u32, f32, f32)] = ~[];
    let mut mint = inmint;
    let mut maxt = inmaxt;
    let mut cur_node = kd_tree_root;
    loop {

        match kd_tree_nodes[cur_node] {
            model::leaf(tri_begin, tri_count) => {

                let mut tri_index = tri_begin;
                while tri_index < tri_begin + tri_count {
                    let t = &get_triangle( polys, tri_index as uint);
                    if ( r.intersect(t).is_some() ) {
                        return true;
                    }
                    tri_index += 1;
                }
                if ( stack.len() > 0 ){
                    let (n,mn,mx) = stack.pop();
                    cur_node = n as uint;
                    mint = mn;
                    maxt = mx;
                } else {
                    return false;
                }
            }
            model::node(axis, splitter, right_tree) => {

                // find the scalar direction/origin for the current axis
                let (inv_dir_scalar, origin) = match axis {
                    model::x() => { (inv_dir.x, r.origin.x) }
                    model::y() => { (inv_dir.y, r.origin.y) }
                    model::z() => { (inv_dir.z, r.origin.z) }
                };

                // figure out which side of the spliting plane the ray origin is
                // i.e. which child we need to test first.
                let (near,far) = if origin < splitter || (origin == splitter && inv_dir_scalar >= 0.0) {
                    ((cur_node+1) as u32,right_tree)
                } else {
                    (right_tree, (cur_node+1) as u32)
                };

                // find intersection with plane
                // origin + dir*t = splitter
                let plane_dist = (splitter - origin) * inv_dir_scalar;

                if plane_dist > maxt || plane_dist < 0.0 {
                    cur_node = near as uint;
                } else if plane_dist < mint {
                    cur_node = far as uint;
                } else{
                    stack.push((far, plane_dist, maxt));
                    cur_node = near as uint;
                    maxt = plane_dist;
                }
            }
        }
    }
}

#[inline]
fn trace_soup( polys: &model::polysoup, r: &Ray) -> Option<(HitResult, uint)>{

    let mut res : Option<(HitResult, uint)> = option::None;

    for uint::range(0, polys.indices.len() / 3) |tri_ix| {
        let tri = &get_triangle( polys, tri_ix);

        let new_hit = r.intersect(tri);

        match (res, new_hit) {
            (option::None(),option::Some(hit)) => {
                res = option::Some((hit, tri_ix));
            }
            (option::Some((old_hit,_)), option::Some(hit))
                        if hit.t < old_hit.t && hit.t > 0.0 => {
                res = option::Some((hit, tri_ix));
            }
            _ => {}
        }
    }

    return res;
}

struct light {
    pos: vec3,
    strength: f32,
    radius: f32,
    color: vec3
}

fn make_light( pos: vec3, strength: f32, radius: f32, color: vec3 ) -> light {
    light{ pos: pos, strength: strength, radius: radius, color: color }
}

#[inline(always)]
fn direct_lighting( lights: &[light], pos: vec3, n: vec3, view_vec: vec3, rnd: &rand_env, depth: uint, occlusion_probe: &fn(vec3) -> bool ) -> vec3 {

    let mut direct_light = vec3(0.0,0.0,0.0);
    for lights.iter().advance |l| {

        // compute shadow contribution
        let mut shadow_contrib = 0.0;
        let num_samples = match depth { 0 => NUM_LIGHT_SAMPLES, _ => 1 };    // do one tap in reflections and GI

        let rot_to_up = rotate_to_up(normalized(sub(pos,l.pos)));
        let shadow_sample_weight = 1.0 / (num_samples as f32);
        do sample_disk(rnd ,num_samples) |u,v| {        // todo: stratify this

            // scale and rotate disk sample, and position it at the light's location
            let sample_pos = add(l.pos,transform(rot_to_up, vec3(u*l.radius,0f32,v*l.radius) ));

            if !occlusion_probe( sub(sample_pos, pos)) {
                shadow_contrib += shadow_sample_weight;
            }
        }

        let light_vec = sub(l.pos, pos);
        let light_contrib =
            if shadow_contrib == 0.0 {
                vec3(0.0, 0.0, 0.0)
            } else {

                let light_vec_n = normalized(light_vec);
                let half_vector = normalized(add(light_vec_n, view_vec));

                let s = dot(n,half_vector);
                let specular = s.pow(&175.0);

                let atten = shadow_contrib*l.strength*(1.0/length_sq(light_vec) + specular*0.05);

                let intensity = atten*dot(n, light_vec_n );

                scale(l.color, intensity)
            };

        direct_light = add(direct_light, light_contrib);
    }

    return direct_light;
}

#[inline]
fn shade(
    pos: vec3, n: vec3, n_face: vec3, r: &Ray, color: vec3, reflectivity: f32, lights: &[light], rnd: &rand_env, depth: uint,
    occlusion_probe: &fn(vec3) -> bool,
    color_probe: &fn(vec3) -> vec3 ) -> vec3 {

    let view_vec = normalized(sub(r.origin, pos));

    // pass in n or n_face for smooth/flat shading
    let shading_normal = if USE_SMOOTH_NORMALS_FOR_DIRECT_LIGHTING { n } else { n_face };

    let direct_light = direct_lighting(lights, pos, shading_normal, view_vec, rnd, depth, occlusion_probe);
    let reflection = sub(scale(shading_normal, dot(view_vec, shading_normal)*2.0), view_vec);
    let rcolor = if reflectivity > 0.001 { color_probe( reflection ) } else { vec3(0.0,0.0,0.0) };

    let mut ambient;
    //ambient = vec3(0.5f32,0.5f32,0.5f32);

    /*let mut ao = 0f32;
    let rot_to_up = rotate_to_up(n_face);
    const NUM_AO_SAMPLES: uint = 5u;
    sample_stratified_2d( rnd, NUM_AO_SAMPLES, NUM_AO_SAMPLES ) { |u,v|
        let sample_vec = transform(rot_to_up, cosine_hemisphere_sample(u,v) );
        //let sample_vec = cosine_hemisphere_sample(u,v);
        if !occlusion_probe( scale(sample_vec, 0.1f32) ) {
                ao += 1f32/((NUM_AO_SAMPLES*NUM_AO_SAMPLES) as f32);
        }
    };
    ambient = scale(ambient,ao); // todo: add ambient color */

    // Final gather GI
    let gi_normal = if USE_SMOOTH_NORMALS_FOR_GI { n } else { n_face };
    ambient = vec3(0.0,0.0,0.0);
    if depth == 0 && NUM_GI_SAMPLES_SQRT > 0 {
        do sample_cosine_hemisphere( rnd, gi_normal ) |sample_vec| {
            ambient = add( ambient, color_probe( sample_vec ) );
        };
        ambient = scale(ambient, 1.0 / (((NUM_GI_SAMPLES_SQRT * NUM_GI_SAMPLES_SQRT) as f32) * f32::consts::pi ));
    }

    lerp( mul(color,add(direct_light, ambient)), rcolor, reflectivity)
}


struct intersection {
    pos: vec3,
    n: vec3,
    n_face: vec3,
    color: vec3,
    reflectivity: f32
}

#[inline]
fn trace_checkerboard( checkerboard_height: f32, r : &Ray, mint: f32, maxt: f32) -> (Option<intersection>, f32) {
    // trace against checkerboard first
    let checker_hit_t = (checkerboard_height - r.origin.y) / r.dir.y;

    // compute checkerboard color, if we hit the floor plane
    if checker_hit_t > mint && checker_hit_t < maxt {

            let pos = add(r.origin,scale(r.dir, checker_hit_t));

            // hacky checkerboard pattern
            let (u,v) = ((pos.x*5.0).floor() as int, (pos.z*5.0).floor() as int);
            let is_white = (u + v) % 2 == 0;
            let color = if is_white { vec3(1.0,1.0,1.0) } else { vec3(1.0,0.5,0.5) };
            let intersection = option::Some( intersection{
                        pos: pos,
                        n: vec3(0.0,1.0,0.0),
                        n_face: vec3(0.0,1.0,0.0),
                        color: color,
                        reflectivity: if is_white {0.3} else {0.0} } );
            (intersection, checker_hit_t)
    } else {
        (option::None, maxt)
    }
}

#[inline]
fn trace_ray( r : &Ray, mesh : &model::mesh, mint: f32, maxt: f32) -> Option<intersection> {

    let use_kd_tree = true;

    let y_size = sub(mesh.bounding_box.max, mesh.bounding_box.min).y;

    // compute checkerboard color, if we hit the floor plane
    let (checker_intersection, new_maxt) = trace_checkerboard(-y_size*0.5,r,mint,maxt);


    // check scene bounding box first
    if !r.aabb_check( new_maxt, mesh.bounding_box ){
        return checker_intersection;
    }

    // trace against scene
    let trace_result = if use_kd_tree {
        trace_kd_tree( &mesh.polys, mesh.kd_tree.nodes, mesh.kd_tree.root, r, recip( r.dir ), mint, new_maxt )
    } else {
        trace_soup( &mesh.polys, r)
    };

    match trace_result {
        option::Some((hit_info, tri_ix)) if hit_info.t > 0.0 => {
            let pos = add( r.origin, scale(r.dir, hit_info.t));

            let (i0,i1,i2) = (  mesh.polys.indices[tri_ix*3  ],
                                mesh.polys.indices[tri_ix*3+1],
                                mesh.polys.indices[tri_ix*3+2] );

            // interpolate vertex normals...
            let n = normalized(
                    add(    scale( mesh.polys.normals[i0], hit_info.barycentric.z),
                    add(    scale( mesh.polys.normals[i1], hit_info.barycentric.x),
                            scale( mesh.polys.normals[i2], hit_info.barycentric.y))));

            // compute face-normal
            let (v0,v1,v2) = (  mesh.polys.vertices[i0],
                                mesh.polys.vertices[i1],
                                mesh.polys.vertices[i2] );
            let n_face = normalized( cross(sub(v1,v0), sub(v2,v0)));

            option::Some( intersection{
                    pos: pos,
                    n: n,
                    n_face: n_face,
                    color: vec3(1.0,1.0,1.0),
                    reflectivity: 0.0 } )
        }
        _ => {
            checker_intersection
        }
    }
}

#[inline]
fn trace_ray_shadow( r: &Ray, mesh: &model::mesh, mint: f32, maxt: f32) -> bool {

    let y_size = sub(mesh.bounding_box.max, mesh.bounding_box.min).y;

    // compute checkerboard color, if we hit the floor plane
    let (checker_intersection, new_maxt) = trace_checkerboard(-y_size*0.5,r,mint,maxt);

    if ( checker_intersection.is_some() ) {
        return true;
    }

    // check scene bounding box first
    if !r.aabb_check( new_maxt, mesh.bounding_box ){
        return false;
    }

    // trace against scene
    trace_kd_tree_shadow( &mesh.polys, mesh.kd_tree.nodes, mesh.kd_tree.root, r, recip( r.dir ), mint, new_maxt )
}


#[inline(always)]
fn get_color( r: &Ray, mesh: &model::mesh, lights: &[light], rnd: &rand_env, tmin: f32, tmax: f32, depth: uint) -> vec3 {
    let theta = dot( vec3(0.0,1.0,0.0), r.dir );
    let default_color = vec3(clamp(1.0-theta*4.0,0.0,0.75)+0.25, clamp(0.5-theta*3.0,0.0,0.75)+0.25, theta);    // fake sky colour

    if depth >= MAX_TRACE_DEPTH {
        return default_color;
    }
    
    match trace_ray( r, mesh, tmin, tmax ) {
        option::Some(intersection{pos,n,n_face,color,reflectivity}) => {
            let surface_origin = add(pos, scale(n_face, 0.000002));

            shade(pos, n, n_face, r, color, reflectivity, lights, rnd, depth,
            |occlusion_vec| {
                let occlusion_ray = &Ray{origin: surface_origin, dir: occlusion_vec};
                trace_ray_shadow(occlusion_ray, mesh, 0.0, 1.0)
            },
            |ray_dir| {
                let reflection_ray = &Ray{origin: surface_origin, dir: normalized(ray_dir)};
                get_color(reflection_ray, mesh, lights, rnd, tmin, tmax, depth + 1)
            })

        }
        _ => { default_color }
    }

}

#[inline]
fn gamma_correct( v : vec3 ) -> vec3 {
    vec3( v.x.pow( &(1.0/2.2) ),
          v.y.pow( &(1.0/2.2) ),
          v.z.pow( &(1.0/2.2) ))
}

struct TracetaskData {
	meshARC: extra::arc::ARC<model::mesh>,
	horizontalFOV: f32,
	width: uint,
    height: uint,
	sample_grid_size: uint,
	height_start: uint,
	height_stop: uint,
	sample_coverage_inv: f32,
	lights: ~[light],
	rnd: ~rand_env
}

#[inline]
fn tracetask(data: ~TracetaskData) -> ~[Color] {
    match data {
        ~TracetaskData {meshARC: meshARC, horizontalFOV: horizontalFOV, width: width,
            height: height, sample_grid_size: sample_grid_size, height_start: height_start,
            height_stop: height_stop, sample_coverage_inv: sample_coverage_inv, lights: lights,
            rnd: rnd} => {
                let mesh = meshARC.get();
            	let mut img_pixels = vec::with_capacity(width);
            	for uint::range( height_start, height_stop ) |row| {
	                for uint::range( 0, width ) |column| {
                        let mut shaded_color = vec3(0.0,0.0,0.0);
                        
                        do sample_stratified_2d(rnd, sample_grid_size, sample_grid_size) |u,v| {
                            let sample = match sample_grid_size { 1 => (0.0,0.0), _ => (u-0.5,v-0.5) };
                            let r = &get_ray(horizontalFOV, width, height, column, row, sample );
                            shaded_color = add( shaded_color, get_color(r, mesh, lights, rnd, 0.0, f32::infinity, 0));
                        }
                        shaded_color = scale(gamma_correct(scale( shaded_color, sample_coverage_inv * sample_coverage_inv)), 255.0);
                        let pixel = Color{  r: clamp(shaded_color.x, 0.0, 255.0) as u8,
                                            g: clamp(shaded_color.y, 0.0, 255.0) as u8,
                                            b: clamp(shaded_color.z, 0.0, 255.0) as u8 };
                        img_pixels.push(pixel)
	                }
                }
            img_pixels
        }
    }
}

pub fn generate_raytraced_image_single(
    mesh: model::mesh,
    horizontalFOV: f32,
    width: uint,
    height: uint,
    sample_grid_size: uint,
    sample_coverage_inv: f32,
    lights: ~[light]) -> ~[Color]
{
    let rnd = get_rand_env();
    for_each_pixel(width, height, |x,y| {
        let mut shaded_color = vec3(0.0,0.0,0.0);
        
        do sample_stratified_2d(&rnd, sample_grid_size, sample_grid_size) |u,v| {
            let sample = match sample_grid_size { 1 => (0.0,0.0), _ => (u-0.5,v-0.5) };
            let r = &get_ray(horizontalFOV, width, height, x, y, sample );
            shaded_color = add( shaded_color, get_color(r, &mesh, lights, &rnd, 0.0, f32::infinity, 0));
        }
        shaded_color = scale(gamma_correct(scale( shaded_color, sample_coverage_inv*sample_coverage_inv)), 255f32);
        Color{
            r: clamp(shaded_color.x, 0.0, 255.0) as u8,
            g: clamp(shaded_color.y, 0.0, 255.0) as u8,
            b: clamp(shaded_color.z, 0.0, 255.0) as u8 }
    })
}

// This fn generates the raytraced image by spawning 'num_tasks' tasks, and letting each
// generate a part of the image. The way the work is divided is not intelligent: it chops
// the image in horizontal chunks of step_size pixels high, and divides these between the
// tasks. There is no work-stealing :(
pub fn generate_raytraced_image_multi(
    mesh: model::mesh,
    horizontalFOV: f32,
    width: uint,
    height: uint,
    sample_grid_size: uint,
    sample_coverage_inv: f32,
    lights: ~[light],
    num_tasks: uint) -> ~[Color]
{
    io::print(fmt!("using %? tasks ... ", num_tasks));
    let meshARC = extra::arc::ARC(mesh);
    let rnd = get_rand_env();
    let mut workers = ~[];
    for num_tasks.times {
        workers.push(concurrent::ConcurrentCalc::new());
    }
    let step_size = 4;
    let mut results = ~[];
    for uint::range(0,(height / step_size)+1) |i| {
        let ttd = ~TracetaskData {   // The data required to trace the rays.
            meshARC: meshARC.clone(),
            horizontalFOV: horizontalFOV,
            width: width,
            height: height,
            sample_grid_size: sample_grid_size,
            height_start: uint::min( i * step_size, height),
            height_stop: uint::min( (i + 1) * step_size, height ),
            sample_coverage_inv: sample_coverage_inv,
            lights: copy lights,
            rnd: ~rnd.clone()
        };
        results.push(workers[i % num_tasks].calculate(ttd,tracetask));
    }
//    println("Done distributing the work, let's evaluate!");
    let mut result = ~[];
    for results.consume_iter().advance |future| {
        let mut future = future;
        result.push_all( future.get() )
    }
    result
}

extern {
    fn rust_num_threads() -> libc::uintptr_t;   // A trick that should tell us the number of processors.
}

pub fn generate_raytraced_image(
    mesh: model::mesh,
    horizontalFOV: f32,
    width: uint,
    height: uint,
    sample_grid_size: uint) -> ~[Color]
{
    let sample_coverage_inv = 1.0 / (sample_grid_size as f32);
    let lights = ~[  make_light(vec3(-3.0, 3.0, 0.0),10.0, 0.3, vec3(1.0,1.0,1.0)) ]; //,
                    //make_light(vec3(0f32, 0f32, 0f32), 10f32, 0.25f32, vec3(1f32,1f32,1.0f32))];
    let mut num_tasks = match NUM_THREADS {
      0 => unsafe { rust_num_threads() as uint },
      n => n
    };
    if num_tasks > height { num_tasks = height };   // We evaluate complete rows, there is no point in having more tasks than there are rows.
    match num_tasks {
        1 => generate_raytraced_image_single(mesh,horizontalFOV,width,height,sample_grid_size,sample_coverage_inv,lights),
        n => generate_raytraced_image_multi(mesh,horizontalFOV,width,height,sample_grid_size,sample_coverage_inv,lights,n)
    }   
}
