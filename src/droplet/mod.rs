//! Droplet-style (10x 3-read) single-cell Hi-C.
mod barcode;
mod config;
mod fragments;
mod index;
mod qc;
mod run;
pub(crate) mod slurm;
pub(crate) mod two_read;

pub(crate) use run::run;
