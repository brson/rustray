import std::{option,math};

type vec3 = { x:float, y:float, z:float};

fn vec3(x:float, y:float, z:float) -> vec3{
	{x:x, y:y, z:z}
}

fn dot(a:vec3, b:vec3) -> float {
	a.x*b.x + a.y*b.y + a.z*b.z
}

fn lerp(a:vec3, b:vec3, t:float) -> vec3{
	add( a, scale(sub(b,a),t) )
}

fn scale(v:vec3, s:float) -> vec3 {
	{ x:v.x*s, y:v.y*s, z:v.z*s }
}

fn length_sq(v:vec3) -> float {
	dot(v,v)
}

fn length(v:vec3) -> float {
	math::sqrt(length_sq(v))
}

fn normalized(v:vec3) -> vec3 {
	scale(v, 1.0f / length(v))
}

fn mul(a:vec3, b:vec3) -> vec3{
	vec3(a.x*b.x, a.y*b.y, a.z*b.z)
}

fn add(a:vec3, b:vec3) -> vec3 {
	{x:a.x+b.x, y:a.y+b.y, z:a.z+b.z}
}

fn sub(a:vec3, b:vec3) -> vec3 {
	add(a, scale(b, -1.0f))
}

fn cross(a:vec3, b:vec3) -> vec3 {
	{ x: a.y*b.z - b.y*a.z, 
	  y: a.z*b.x - b.z*a.x, 
	  z: a.x*b.y - b.x*a.y }
}

fn min(a: vec3, b: vec3) -> vec3 {
	vec3( math::min(a.x,b.x), math::min(a.y, b.y), math::min(a.z, b.z) )
}

fn max(a: vec3, b: vec3) -> vec3 {
	vec3( math::max(a.x,b.x), math::max(a.y, b.y), math::max(a.z, b.z) )
}

type ray = { origin:vec3, dir:vec3 };

type triangle = { p1: vec3, p2: vec3, p3: vec3 };

type hit_result = { barycentric: vec3, t: float };

fn ray_triangle_intersect( r:ray, t:triangle ) -> option::t<hit_result> {
	let e1 = sub(t.p2, t.p1);
	let e2 = sub(t.p3, t.p1);
	let s1 = cross(r.dir, e2);
	let divisor = dot(s1,e1);

	if divisor == 0f {
		ret option::none;
	}

	// compute first barycentric coordinate
	let inv_divisor = 1.0f / divisor;
	let d = sub(r.origin, t.p1);
	
	let b1 = dot(d, s1) * inv_divisor;
	if b1 < 0.0f || b1 > 1.0f {
		ret option::none; // outside triangle
	}

	// and second barycentric coordinate
	let s2 = cross(d, e1);
	let b2 = dot(r.dir, s2)*inv_divisor;
	if b2 < 0.0f || b1+b2 > 1.0f {
		ret option::none; // outside triangle
	}

	let t = dot(e2, s2) * inv_divisor;
	if t < 0.0f {
		ret option::none; // behind viewer
	}
	
	option::some( {barycentric: { x:b1, y: b2, z: 1.0f-b1-b2 }, t: t} )
}

#[test]
fn intersection_test()
{
	let ray : ray = { origin: vec3(0f, 0f, 0f), dir: vec3(0.0f,0.0f,-1.0f) };
	let tri : triangle = { 	p1: vec3(-1f, -1f, -1f), 
					p2: vec3(1f, -1f, -1f), 
					p3: vec3(0f, 2f, -1f) };

	assert option::is_some( ray_triangle_intersect( ray, tri ) );
}
