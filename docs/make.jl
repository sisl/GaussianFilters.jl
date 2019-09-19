using Documenter, GaussianFilters

# This function builds the documentation
makedocs(
    modules   = [GaussianFilters],
    doctest   = false,
    #clean     = true,
    #linkcheck = false,
    format    = Documenter.HTML(),
    sitename  = "GaussianFilters.jl",
    pages     =
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "Introduction" => [
            "Basics" => "index.md",
            "Installation" => "install.md",
           ],

        "User Documentation" => [
            "Kalman-class Filters" => "kalman.md",
            "GM-PHD Filter" => "gmphd.md",
           ]
    ]
)

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
