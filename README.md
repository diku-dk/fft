# Fast Fourier Transforms in Futhark [![Build Status](https://travis-ci.org/diku-dk/fft.svg?branch=master)](https://travis-ci.org/diku-dk/fft) [![Documentation](https://futhark-lang.org/pkgs/github.com/athas/fft/status.svg)](https://futhark-lang.org/pkgs/github.com/athas/fft/latest/)

A library for perfoming FFTs in Futhark.  Currently only provides a
radix-2 Stockham implementation, but will hopefully grow in the future
to contain more algorithms.

## Installation

```
$ futhark pkg add github.com/diku-dk/fft
$ futhark pkg sync
```

## Usage

```
> import "lib/github.com/diku-dk/fft/stockham-radix-2"
> module fft32 = mk_fft f32
> unzip (fft32.fft (zip [1,2,3] [4,5,6]))
[6.0f32, 3.0f32, 2.0f32]
[15.0f32, -4.0f32, 5.0f32]
```

## See also

* https://github.com/diku-dk/complex
