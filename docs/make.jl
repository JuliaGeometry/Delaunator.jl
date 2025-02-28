using Documenter, Delaunator

makedocs(;
    modules=[Delaunator],
    format=Documenter.HTML(; description = "A port of the Delaunator algorithm to Julia"),
    pages=[
        "Home" => "index.md",
        "API" => "api.md"
    ],
    repo="https://github.com/JuliaGeometry/Delaunator.jl/blob/{commit}{path}#L{line}",
    sitename="Delaunator.jl",
    authors="David Gleich, Steve Kelly",
)

deploydocs(;
    repo="github.com/JuliaGeometry/Delaunator.jl",
)
