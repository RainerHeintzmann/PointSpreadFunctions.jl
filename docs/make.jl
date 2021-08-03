using Base: package_slug
using Documenter, PSFs

makedocs(modules=[PSFs], sitename="PSFs Documentation", doctest = false,
        pages = Any["PSFs.jl" => "index.md",
        "Workflow" => Any[
                        "workflow/PSF_parameters.md",
                        "workflow/PSF_generation.md",
                    ],
        "Function references" => Any[
                    "function_references/all_functions.md",
                    ]
                    ]
        )

deploydocs(repo = "github.com/RainerHeintzmann/PSFs.jl.git", devbranch = "main")
