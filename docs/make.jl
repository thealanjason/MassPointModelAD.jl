using MPMad
using Documenter

DocMeta.setdocmeta!(MPMad, :DocTestSetup, :(using MPMad); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [MPMad],
    authors = "Alan Correa <correa@mbd.rwth-aachen.de>",
    repo = "https://github.com/thealanjason/MPMad.jl/blob/{commit}{path}#{line}",
    sitename = "MPMad.jl",
    format = Documenter.HTML(; canonical = "https://thealanjason.github.io/MPMad.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/thealanjason/MPMad.jl")
