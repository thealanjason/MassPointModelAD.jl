using MassPointModelAD
using Documenter

DocMeta.setdocmeta!(MassPointModelAD, :DocTestSetup, :(using MassPointModelAD); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [MassPointModelAD],
    authors = "Alan Correa <correa@mbd.rwth-aachen.de>",
    repo = "https://github.com/thealanjason/MassPointModelAD.jl/blob/{commit}{path}#{line}",
    sitename = "MassPointModelAD.jl",
    format = Documenter.HTML(; canonical = "https://thealanjason.github.io/MassPointModelAD.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/thealanjason/MassPointModelAD.jl")
