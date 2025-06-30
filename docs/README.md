# NonlinearNormalForm Documentation Cheat Sheet

## Documenter

The NonlinearNormalForm documentation is built using Documenter. Documenter documentation at:
- [https://documenter.juliadocs.org/stable/man/guide/](https://documenter.juliadocs.org/stable/man/guide/)

Some notes on how to use Documenter at:
- [https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/](https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/)


## Building the Documentation Locally

Run the following command from the docs/ directory
```
$ julia --project make.jl
```
Note that $ just represents the prompt character. You don't need to type that.
The output you see should be something like:
```
(base) MAC:~/.julia/dev/NonlinearNormalForm/docs>  julia --project make.jl
[ Info: SetupBuildDirectory: setting up build directory.
[ Info: Doctest: running doctests.
[ Info: ExpandTemplates: expanding markdown templates.
[ Info: CrossReferences: building cross-references.
[ Info: CheckDocument: running document checks.
[ Info: Populate: populating indices.
[ Info: RenderDocument: rendering document.
[ Info: HTMLWriter: rendering HTML pages.
[ Info: Automatic `version="0.2.0"` for inventory from ../Project.toml
┌ Warning: Documenter could not auto-detect the building environment. Skipping deployment.
└ @ Documenter ~/.julia/packages/Documenter/uNI6f/src/deployconfig.jl:93
```
The warning is fine since you are just building the docs locally.

## Visualize the Docs locally:

Do the following:
```
julia> ] activate docs    # Do this in the NonlinearNormalForm directory
julia> using LiveServer 
julia> servedocs()
```
and the docs will be rendered and hosted locally at the URL provided in the output.
Generally the URL will be something like [http://localhost:8000/](http://localhost:8000/)
