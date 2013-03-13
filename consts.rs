// TODO: These things should all be overridable by the command line
pub const WIDTH : uint = 1024u;
pub const HEIGHT : uint = 1024u;
pub const FOV : f32 = 3.14159f32 / 3f32 ;
pub const SAMPLE_GRID_SIZE : uint = 1u;
pub const NUM_GI_SAMPLES_SQRT: uint = 4u;
pub const NUM_LIGHT_SAMPLES : uint = 8u;
pub const MAX_TRACE_DEPTH : uint = 1u;
pub const USE_SMOOTH_NORMALS_FOR_GI : bool = true;
pub const USE_SMOOTH_NORMALS_FOR_DIRECT_LIGHTING : bool = true;
pub const NUM_THREADS: uint = 32u;
