# NonlinearNormalForm Documentation Cheat Sheet

## Documenter

The NonlinearNormalForm documentation is built using Documenter. Documenter documentation at:
- [https://documenter.juliadocs.org/stable/man/guide/](https://documenter.juliadocs.org/stable/man/guide/)

Some notes on how to use Documenter at:
- [https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/](https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/)


## Building the Documentation Locally

First download `NonlinearNormalForm` to the `~/.julia/dev/` area. Then make sure that you are
running julia with `NonlinearNormalForm` in `development` mode 
[https://pkgdocs.julialang.org/v1/managing-packages/#Adding-a-local-package](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-a-local-package).
This is important since you want to have local changes reflected in the documentation build so
that you can visualize your changes locally.

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

Note: You do not need to build the documentation if you use the server (see the next section)
as the server will automatically do this.

## Visualize the Docs locally:

Do the following in the `NonlinearNormalForm` directory:
```
julia> ] activate docs
julia> using LiveServer 
julia> servedocs()
```
and the docs will be rendered and hosted locally at the URL provided in the output.
Generally the URL will be something like [http://localhost:8000/](http://localhost:8000/)

Note: If you use the server, the server will build documentation as needed. That is, the documentation
is rebuilt whenever there are changes to the markdown files. Also, when using the server, you
do not have to build the docs as outlined in the last section.

If you suspect there are detached server processes that are not closing, use the following command to list them:
```
$ lsof | grep julia
```
And then use `kill -9 <process-number>` to kill them where `<process-number>` is given in the
second column of the `lsof` output.