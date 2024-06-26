% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_surrogate.R
\name{rand_fft}
\alias{rand_fft}
\title{func_surrogate}
\usage{
rand_fft(amp, n, sigma)
}
\arguments{
\item{amp}{Numeric, the amplitude of target time series.}

\item{n}{Integer, the length of target time series.}

\item{sigma}{Numeric, the standard deviation of target time series.}
}
\description{
R function file for 'andsr' (Analysis of Nonlinear Dynamical Systems in R).
This file contains functions related to surrogate data generation.
Intially written on 20200227 by K.Kawatsu.
Last update: 20220502.
FFT randomization

\code{rand_fft} returns a time-series vector shuffled with
FFT (Fast-Fourier Transformation) randomization.
}
