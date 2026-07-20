using Documenter, GaussianFilters

makedocs(
    modules   = [GaussianFilters],
    format    = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true"),
    sitename  = "GaussianFilters",
    pages     = ["Introduction" => [
                    "Basics" => "index.md",
                    "Installation" => "install.md"
                    ],
                "User Documentation" => [
                    "Kalman-class Filters" => "kalman.md",
                    "GM-PHD Filter" => "gmphd.md"
                    ],
                "Integrations" => [
                    "POMDPs.jl" => "pomdps.md"
                    ]
                ])

if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo = "github.com/sisl/GaussianFilters.jl",
        push_preview = true,
    )
end
