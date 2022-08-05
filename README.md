# PointSpreadFunctions.jl
Toolbox for calculating optical PointSpreadFunctions (PSFs). It incorporates aberrations (specified though Zernike coefficients) and supports a number of
models for calculation (`PointSpreadFunctions.MethodSincR`, `PointSpreadFunctions.MethodPropagate`, `PointSpreadFunctions.MethodPropagateIterative`), some are faster and some more accurate.

For detailed documentation, see the docs/dev (or docs/stable) link below.
You may also find the file `PSFs_demo.ipnb` in the examples folder helpful.

## Installation
To get the latest stable release of PointSpreadFunctions.jl type `]` in the Julia REPL:
```
] add PointSpreadFunctions.jl
```


| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |



## JuliaCon 2022
<a  href="https://www.youtube.com/watch?v=HwhrCIGzYAk"><img src="docs/src/assets/juliacon.png"  width="300"></a>




[docs-dev-img]: https://img.shields.io/badge/docs-dev-orange.svg 
[docs-dev-url]: https://RainerHeintzmann.github.io/PointSpreadFunctions.jl/dev/ 

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg 
[docs-stable-url]: https://RainerHeintzmann.github.io/PointSpreadFunctions.jl/stable/

[codecov-img]: https://codecov.io/gh/RainerHeintzmann/PointSpreadFunctions.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/RainerHeintzmann/PointSpreadFunctions.jl

[CI-img]: https://github.com/RainerHeintzmann/PointSpreadFunctions.jl/workflows/CI/badge.svg
[CI-url]: https://github.com/RainerHeintzmann/PointSpreadFunctions.jl/actions?query=workflow%3ACI 
