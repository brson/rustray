

type vec3 = { x:f32, y:f32, z:f32 };
type mtx33 = { r0:vec3, r1:vec3, r2:vec3 };

#[inline(always)] 
fn vec3(x:f32, y:f32, z:f32) -> vec3{
	{x:x, y:y, z:z}
}

#[inline(always)] 
fn dot(a:vec3, b:vec3) -> f32 {
	a.x*b.x + a.y*b.y + a.z*b.z
}

#[inline(always)] 
fn lerp(a:vec3, b:vec3, t:f32) -> vec3{
	add( a, scale(sub(b,a),t) )
}

#[inline(always)] 
fn scale(v:vec3, s:f32) -> vec3 {
	{ x:v.x*s, y:v.y*s, z:v.z*s }
}

#[inline(always)] 
fn length_sq(v:vec3) -> f32 {
	dot(v,v)
}

#[inline(always)] 
fn length(v:vec3) -> f32 {
	f32::sqrt(length_sq(v))
}

#[inline(always)] 
fn normalized(v:vec3) -> vec3 {
	scale(v, 1.0f32 / length(v))
}

#[inline(always)] 
fn recip(a:vec3) -> vec3{
	vec3(1f32/a.x, 1f32/a.y, 1f32/a.z)
}


#[inline(always)] 
fn mul(a:vec3, b:vec3) -> vec3{
	vec3(a.x*b.x, a.y*b.y, a.z*b.z)
}

#[inline(always)] 
fn add(a:vec3, b:vec3) -> vec3 {
	{x:a.x+b.x, y:a.y+b.y, z:a.z+b.z}
}

#[inline(always)] 
fn sub(a:vec3, b:vec3) -> vec3 {
	add(a, scale(b, -1.0f32))
}

#[inline(always)] 
fn cross(a:vec3, b:vec3) -> vec3 {
	{ x: a.y*b.z - b.y*a.z, 
	  y: a.z*b.x - b.z*a.x, 
	  z: a.x*b.y - b.x*a.y }
}

#[inline(always)] 
fn min(a: vec3, b: vec3) -> vec3 {
	vec3( f32::fmin(a.x,b.x), f32::fmin(a.y, b.y), f32::fmin(a.z, b.z) )
}

#[inline(always)] 
fn max(a: vec3, b: vec3) -> vec3 {
	vec3( f32::fmax(a.x,b.x), f32::fmax(a.y, b.y), f32::fmax(a.z, b.z) )
}

type ray = { origin:vec3, dir:vec3 };

type triangle = { p1: vec3, p2: vec3, p3: vec3 };

type hit_result = { barycentric: vec3, t: f32 };

#[inline(always)] 
fn ray_triangle_intersect( r:ray, t:triangle ) -> option<hit_result> {
	let e1 = sub(t.p2, t.p1);
	let e2 = sub(t.p3, t.p1);
	let s1 = cross(r.dir, e2);
	let divisor = dot(s1,e1);

	if divisor == 0f32 {
		ret option::none;
	}

	// compute first barycentric coordinate
	let inv_divisor = 1.0f32 / divisor;
	let d = sub(r.origin, t.p1);
	
	let b1 = dot(d, s1) * inv_divisor;
	if b1 < 0.0f32 || b1 > 1.0f32 {
		ret option::none; // outside triangle
	}

	// and second barycentric coordinate
	let s2 = cross(d, e1);
	let b2 = dot(r.dir, s2)*inv_divisor;
	if b2 < 0.0f32 || b1+b2 > 1.0f32 {
		ret option::none; // outside triangle
	}

	let t = dot(e2, s2) * inv_divisor;
	if t < 0.0f32 {
		ret option::none; // behind viewer
	}
	
	option::some( {barycentric: { x:b1, y: b2, z: 1.0f32-b1-b2 }, t: t} )
}

type aabb = { min: vec3, max: vec3 };

#[inline(always)] 
fn ray_aabb_check( r:ray, max_dist: f32, box: aabb ) -> bool {
	let inv_dir = recip(r.dir);
	let (tx1,tx2,ty1,ty2,tz1,tz2) = (
		(box.min.x - r.origin.x)*inv_dir.x,
		(box.max.x - r.origin.x)*inv_dir.x,
		(box.min.y - r.origin.y)*inv_dir.y,
		(box.max.y - r.origin.y)*inv_dir.y,
		(box.min.z - r.origin.z)*inv_dir.z,
		(box.max.z - r.origin.z)*inv_dir.z
	);
	
	let (minx, maxx) = (f32::fmin(tx1,tx2), f32::fmax(tx1,tx2));
	let (miny, maxy) = (f32::fmin(ty1,ty2), f32::fmax(ty1,ty2));
	let (minz, maxz) = (f32::fmin(tz1,tz2), f32::fmax(tz1,tz2));
	
	let tmin = f32::fmax(minx, f32::fmax(miny, minz));
	let tmax = f32::fmin(maxx, f32::fmin(maxy, maxz));

	tmax >= 0f32 && tmin <= tmax && tmin <= max_dist
}

// Gives a cosine hemisphere sample from two uniform f32s
// in [0,1) range.
#[inline(always)] 
fn cosine_hemisphere_sample( u: f32, v: f32 ) -> vec3 {
	let r_sqrt = f32::sqrt(u);
	let theta = 2f32 * f32::consts::pi * v;
	vec3( r_sqrt*f32::cos(theta), f32::sqrt(1f32-u), r_sqrt*f32::sin(theta))
}

#[inline(always)] 
fn rotate_to_up( up_vec: vec3 ) -> mtx33 {
	let perp = if up_vec == vec3(0f32,1f32,0f32) { vec3(1f32,0f32,0f32) } else { vec3(0f32,1f32,0f32) };
	let right = cross( up_vec, perp );
	let fwd = cross(right, up_vec );	
	transposed( { r0: right, r1: up_vec, r2: fwd } )
}

#[inline(always)]
fn rotate_y(theta: f32) -> mtx33{
	let ct = f32::cos(theta);
	let st = f32::sin(theta);
	{ r0: vec3(ct,0f32,st), r1: vec3(0f32,1f32,0f32), r2: vec3(-st, 0f32, ct) }
}

#[inline(always)] 
fn transform( m: mtx33, v: vec3 ) -> vec3 {
	vec3( dot( m.r0, v ), dot( m.r1, v ), dot( m.r2, v ) )
}

#[inline(always)] 
fn mul_mtx33( a: mtx33, b: mtx33 ) -> mtx33 {
	let b_t = transposed(b);
	{ 	r0: vec3( dot(a.r0,b_t.r0), dot(a.r0,b_t.r1), dot(a.r0,b_t.r2) ),
	 	r1: vec3( dot(a.r1,b_t.r0), dot(a.r1,b_t.r1), dot(a.r1,b_t.r2) ),
	 	r2: vec3( dot(a.r2,b_t.r0), dot(a.r2,b_t.r1), dot(a.r2,b_t.r2) ) }
}

#[inline(always)] 
fn transposed( m: mtx33 ) -> mtx33 {
	{ 	r0: vec3( m.r0.x, m.r1.x, m.r2.x ),
		r1: vec3( m.r0.y, m.r1.y, m.r2.y ),
		r2: vec3( m.r0.z, m.r1.z, m.r2.z ) }
}

#[test]
fn intersection_test()
{
	let ray : ray = { origin: vec3(0f32, 0f32, 0f32), dir: vec3(0.0f32,0.0f32,-1.0f32) };
	let tri : triangle = { 	p1: vec3(-1f32, -1f32, -1f32), 
					p2: vec3(1f32, -1f32, -1f32), 
					p3: vec3(0f32, 2f32, -1f32) };

	assert option::is_some( ray_triangle_intersect( ray, tri ) );
}
