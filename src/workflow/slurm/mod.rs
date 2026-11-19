//! Slurm job-array execution: index, parallel alignment, and final processing.
mod state;
mod status;
mod submit;
mod worker;

pub(crate) use state::{clear_failure, now_ms, record_failure, write_json, Jobs, SCHEMA};
pub(crate) use status::{cancel, status};
pub(crate) use submit::{
    cancel_submitted, release_array, sbatch_named, submit, validate_extra_args,
    write_task_script,
};
pub(crate) use worker::task;
