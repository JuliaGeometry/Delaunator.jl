using Documenter, Delaunator

makedocs(;
    modules=[Delaunator],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/sjkelly/Delaunator.jl/blob/{commit}{path}#L{line}",
    sitename="Delaunator.jl",
    authors="Steve Kelly",
    assets=String[],
)

deploydocs(;
    repo="github.com/sjkelly/Delaunator.jl",
)
