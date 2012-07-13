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
