// TODO: These things should all be overridable by the command line
const WIDTH : uint = 128u;
const HEIGHT : uint = 128u;
const FOV : f32 = 3.14159f32 / 3f32 ;
const SAMPLE_GRID_SIZE : uint = 1u;
const NUM_GI_SAMPLES_SQRT: uint = 4u;
const NUM_LIGHT_SAMPLES : uint = 8u;
const MAX_TRACE_DEPTH : uint = 1u;
const USE_SMOOTH_NORMALS_FOR_GI : bool = true;
const USE_SMOOTH_NORMALS_FOR_DIRECT_LIGHTING : bool = true;
