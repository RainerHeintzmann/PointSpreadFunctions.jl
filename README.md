# PSFs.jl
Toolbox for calculating optical PSFs. It incorporates aberrations (specified though Zernike coefficients) and supports a number of
models for calculation (`PSFs.MethodSincR`, `PSFs.MethodPropagate`, `PSFs.MethodPropagateIterative`), some are faster and some more accurate.

For detailed documentation, see the docs/dev (or docs/stable) link below.
You may also find the file `PSFs_demo.ipnb` in the examples folder helpful.

## Installation
To get the latest stable release of PSFs.jl type `]` in the Julia REPL:
```
] add https://github.com/RainerHeintzmann/PSFs.jl
```


| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |



[docs-dev-img]: https://img.shields.io/badge/docs-dev-orange.svg 
[docs-dev-url]: https://RainerHeintzmann.github.io/PSFs.jl/dev/ 

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg 
[docs-stable-url]: https://RainerHeintzmann.github.io/PSFs.jl/stable/

[codecov-img]: https://codecov.io/gh/RainerHeintzmann/PSFs.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/RainerHeintzmann/PSFs.jl

[CI-img]: https://github.com/RainerHeintzmann/PSFs.jl/workflows/CI/badge.svg
[CI-url]: https://github.com/RainerHeintzmann/PSFs.jl/actions?query=workflow%3ACI 
