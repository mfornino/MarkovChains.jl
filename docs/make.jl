using MarkovChains
using Documenter

makedocs(;
    modules=[MarkovChains],
    authors="Michele Fornino <mfornino@mit.edu>",
    repo="https://github.com/mfornino/MarkovChains.jl/blob/{commit}{path}#L{line}",
    sitename="MarkovChains.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfornino.github.io/MarkovChains.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="mfornino.github.io/MarkovChains.jl/",
)
