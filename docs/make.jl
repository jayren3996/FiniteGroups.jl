using Documenter
using FiniteGroups
using LinearAlgebra

DocMeta.setdocmeta!(FiniteGroups, :DocTestSetup, :(using FiniteGroups, LinearAlgebra); recursive=true)

makedocs(;
    sitename="FiniteGroups.jl",
    modules=[FiniteGroups],
    checkdocs=:exports,
    authors="jayren3996",
    format=Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true"),
    pages=[
        "Overview" => "index.md",
        "Guide" => [
            "Getting Started" => "getting-started.md",
            "Constructing Groups" => "groups.md",
            "Character Tables" => "character-tables.md",
            "Representations" => "representations.md",
        ],
        "Reference" => [
            "API Reference" => "api.md",
        ],
    ],
)

deploydocs(
    repo="github.com/jayren3996/FiniteGroups.jl.git",
    devbranch="main",
)
