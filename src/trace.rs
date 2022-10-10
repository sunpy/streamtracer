//! Streamline tracing functionality.
use num_derive::ToPrimitive;
use numpy::ndarray::{Array, Array1, Array3, ArrayView1, ArrayView2, ArrayView4, ArrayViewMut2, s};

use crate::field::{VectorField, Bounds};

/// Enum denoting status of the streamline tracer
#[derive(PartialEq, Debug, ToPrimitive)]
pub enum TracerStatus{
    /// Still running
    Running,
    /// Ran out of steps
    RanOutOfSteps = 1,
    /// Stepped out of bounds
    OutOfBounds = 2,
}

/// A single stream line status
pub struct StreamlineStatus{
    /// Reason of termination
    pub rot: TracerStatus,
    /// Number of valid coordinates returned
    /// Can be used to slice away extra bits of the array that were not
    /// used for tracing using `xs.slice(s![:n_points, ..])`.
    pub n_points: usize
}

/// Trace streamlines
///
/// # Parameters
///
/// * `seeds` - Seed points for streamlines. Must be shape (nseeds, 3).
/// * `field` - Vector field to track through.
/// * `direction` - Direction to trace in, `1` for forwards, `-1` for backwards.
/// * `step_size` - Size of each individual step to take.
/// * `max_steps` - Maximum number of steps to take per streamline.
///   This directly sets the size of the output streamline array.
pub fn trace_streamlines<'a>(
    seeds: ArrayView2<'a, f64>,
    xgrid: ArrayView1<'a, f64>,
    ygrid: ArrayView1<'a, f64>,
    zgrid: ArrayView1<'a, f64>,
    values: ArrayView4<'a, f64>,
    cyclic: ArrayView1<'a, bool>,
    direction: i32,
    step_size: f64,
    max_steps: usize,
) ->  (Vec<StreamlineStatus>, Array3<f64>) {
    let field = VectorField::new(xgrid, ygrid, zgrid, values, cyclic);
    let n_seeds: usize = seeds.shape()[0];

    let mut xs = Array::zeros((n_seeds, max_steps, 3));
    let mut statuses: Vec<StreamlineStatus> = vec![];

    // Trace from each seed in turn
    for i in 0..n_seeds{
        statuses.push(
            trace_streamline(
                seeds.slice(s![i, ..]),
                &field,
                &direction,
                &step_size,
                xs.slice_mut(s![i, .., ..])
            )
        )
    }

    return (statuses, xs)
}

/// Trace a single streamline
///
/// # Parameters
/// - `x0`: Streamline seed point.
/// - `field`: Vector field to trace through
/// - `direction`: Direction to trace in. Can be 1 for forwards or -1 for backwards.
/// - `step_size`: Step size to take.
/// - `xs`: Output array. Maximum number of steps to take is set by this.
pub fn trace_streamline(
    x0: ArrayView1<f64>,
    field: &VectorField,
    direction: &i32,
    step_size: &f64,
    mut xs: ArrayViewMut2<f64>,
) -> StreamlineStatus{
    // Tracer status
    let mut status = TracerStatus::Running;
    // Number of points traced
    let mut n_points: usize = 1;
    // Fold direction into the definition of step size,
    // using sign(step_size) to determine step direction
    let step = (*step_size) * (*direction as f64);
    let max_steps: usize = xs.shape()[0];

    // Create output array
    // Take a copy of input seed
    let mut x = Array::zeros(3);
    x.assign(&x0);

    // Take streamline steps
    for i in 0..max_steps{
        // Copy current coordinate
        for j in 0..3 {
            xs[[i, j]] = x[[j]];
        }
        // +1 to account for the initial point
        n_points = i+1;
        // Take a single step
        // Updates `x` in place.
        x = rk4_update(x, field, &step);
        x = field.wrap_cyclic(x);
        // Check new point isn't out of bounds
        match field.check_bounds(x.view()){
            Bounds::Out => {
                status = TracerStatus::OutOfBounds;
                break;
            }
            Bounds::In => {
                continue;
            }
        }
    }

    // Finished tracing maximum steps
    if status == TracerStatus::Running {
        status = TracerStatus::RanOutOfSteps;
    }

    return StreamlineStatus{
        rot: status,
        n_points
    };
}

// Update a coordiante (`x`) by taking a single RK4 step
fn rk4_update(
    mut x: Array1<f64>,
    field: &VectorField,
    step_size: &f64
) -> Array1<f64> {
    let mut xu = x.clone();
    let k1 = stream_function(xu.view(), field, step_size);

    xu = x.clone() + 0.5 * k1.clone();
    let k2 = stream_function(xu.view(), field, step_size);

    xu = x.clone() + 0.5 * k2.clone();
    let k3 = stream_function(xu.view(), field, step_size);

    xu = x.clone() + k3.clone();
    let k4 = stream_function(xu.view(), field, step_size);

    let step = (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
    x = x + step;
    return x;
}

/// Return the step that a linear tracing method would take
/// at a given position.
fn stream_function(
    x: ArrayView1<f64>,
    field: &VectorField,
    step_size: &f64
) -> Array1<f64> {
    let vec = field.vector_at_position(x);
    let vmag = (
        vec[[0]].powf(2.) +
        vec[[1]].powf(2.) +
        vec[[2]].powf(2.)
    ).sqrt();
    return (*step_size) * vec / vmag;
}
