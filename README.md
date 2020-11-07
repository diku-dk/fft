# Fast Fourier Transforms in Futhark [![CI](https://github.com/diku-dk/fft/workflows/CI/badge.svg)](https://github.com/diku-dk/fft/actions) [![Documentation](https://futhark-lang.org/pkgs/github.com/diku-dk/fft/status.svg)](https://futhark-lang.org/pkgs/github.com/diku-dk/fft/latest/)

A library for perfoming FFTs in Futhark.  Currently only provides a
radix-2 Stockham implementation, but will hopefully grow in the future
to contain more algorithms.

**Warning:** currently supports *only* input sizes that are powers of
two.

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
