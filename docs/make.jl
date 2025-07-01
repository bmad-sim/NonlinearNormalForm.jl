using NonlinearNormalForm
using Documenter

DocMeta.setdocmeta!(NonlinearNormalForm, :DocTestSetup, :(using NonlinearNormalForm); recursive=true)

makedocs(;
    modules=[NonlinearNormalForm],
    authors="Matt Signorelli",
    sitename="NonlinearNormalForm.jl",
    format=Documenter.HTML(;
        canonical="https://bmad-sim.github.io/NonlinearNormalForm.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => ["index.md", "namespace.md"],
    ],
)

deploydocs(;
    repo="github.com/bmad-sim/NonlinearNormalForm.jl",
    devbranch="main",
)
