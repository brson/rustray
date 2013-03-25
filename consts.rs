// TODO: These things should all be overridable by the command line
pub static WIDTH : uint = 1024u;
pub static HEIGHT : uint = 1024u;
pub static FOV : f32 = 3.14159f32 / 3f32 ;
pub static SAMPLE_GRID_SIZE : uint = 1u;
pub static NUM_GI_SAMPLES_SQRT: uint = 4u;
pub static NUM_LIGHT_SAMPLES : uint = 8u;
pub static MAX_TRACE_DEPTH : uint = 1u;
pub static USE_SMOOTH_NORMALS_FOR_GI : bool = true;
pub static USE_SMOOTH_NORMALS_FOR_DIRECT_LIGHTING : bool = true;
pub static NUM_THREADS: uint = 32u;
