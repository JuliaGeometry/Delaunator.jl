using Documenter, Delaunator

makedocs(;
    modules=[Delaunator],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/juliageometry/Delaunator.jl/blob/{commit}{path}#L{line}",
    sitename="Delaunator.jl",
    authors="David Gleich, Steve Kelly",
    assets=String[],
)

deploydocs(;
    repo="github.com/juliageometry/Delaunator.jl",
)
