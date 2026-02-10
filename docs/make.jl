using LienardWiechert
using Documenter

DocMeta.setdocmeta!(LienardWiechert, :DocTestSetup, :(using LienardWiechert); recursive=true)

makedocs(;
    modules=[LienardWiechert],
    authors="Wolfgang Hogger",
    sitename="LienardWiechert.jl",
    format=Documenter.HTML(;
        canonical="https://howbgl.github.io/LienardWiechert.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/howbgl/LienardWiechert.jl",
    devbranch="main",
)
