using Documenter, GaussianFilters

# This function builds the documentation
makedocs(
    modules   = [GaussianFilters],
    format    = Documenter.HTML(),
    sitename  = "GaussianFilters",
    warnonly  = [:docs_block, :missing_docs],
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
