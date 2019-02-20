using Documenter, JuliaPackageTemplate

#
include("src/makeplots.jl")

# This function builds the documentation
makedocs(
    modules   = [JuliaPackageTemplate],
    doctest   = false,
    clean     = true,
    linkcheck = false,
    format    = Documenter.HTML(),
    sitename  = "JuliaPackageTemplate.jl",
    authors   = "Duncan Eddy",
    pages     = Any[
        "Home" => "index.md",
        "Modules" => Any[
            "modules/submodule.md", # Use default module name in sidebar
            "Rubber Ducks" => "modules/rubber_ducks.md", # Rename a module in sidebar
        ],
        "Examples" => Any[
            "Plotting Example" => "examples/plotting_example.md"
        ],
        "Library Index" => "library_index.md",
    ]
)

# Generate plots
# Note: Must be called after makedocs so the build folder are created
makeplots()

deploydocs(
    repo = "github.com/sisl/JuliaPackageTemplate.jl",
    devbranch = "master",
    devurl = "latest",
    deps = makeplots,
)