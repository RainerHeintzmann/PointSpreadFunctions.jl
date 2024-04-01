using Base: package_slug
using Documenter, PointSpreadFunctions

makedocs(modules=[PointSpreadFunctions], sitename="PSFs Documentation", doctest = false,
        pages = Any["PointSpreadFunctions.jl" => "index.md",
        "Workflow" => Any[
                        "workflow/PSF_parameters.md",
                        "workflow/PSF_generation.md",
                    ],
        "Function references" => Any[
                    "function_references/all_functions.md",
                    ]
                    ],
        warnonly=true,
        )

deploydocs(repo = "github.com/RainerHeintzmann/PointSpreadFunctions.jl.git", devbranch = "main")
