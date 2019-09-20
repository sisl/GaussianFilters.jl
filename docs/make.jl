using Documenter, GaussianFilters

# This function builds the documentation
makedocs(
    modules   = [GaussianFilters],
    #doctest   = false,
    #clean     = true,
    #linkcheck = false,
    format    = :html,
    sitename  = "GaussianFilters",
    pages     = ["Introduction" => [
                    "Basics" => "index.md",
                    "Installation" => "install.md"
                    ],
                "User Documentation" => [
                    "Kalman-class Filters" => "kalman.md",
                    "GM-PHD Filter" => "gmphd.md"
                    ]
                ])

# Generate plots
# Note: Must be called after makedocs so the build folder are created
# makeplots()

#=
deploydocs(
    repo = "github.com/sisl/GaussianFilters.jl",
    devbranch = "master",
    devurl = "latest",
    deps = ,
)
=#
