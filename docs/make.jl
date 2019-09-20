using Documenter, GaussianFilters

# This function builds the documentation
makedocs(
    modules   = [GaussianFilters],
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

deploydocs(
    repo = "github.com/sisl/GaussianFilters.jl",
)
