using Documenter
using DocumenterLaTeX

using LandauGinzburgCategories
using LandauGinzburgCategories.Library

config = (
    modules  = [
        LandauGinzburgCategories,
        LandauGinzburgCategories.Library,
    ],
    repo     = "https://github.com/tkluck/LandauGinzburgCategories.jl.git",
    doctest  = true,
    sitename = "LandauGinzburgCategories.jl",
    authors  = "Timo Kluck",
    pages    = [
        # keep in sync with index.md
        "Home"                => "index.md",
        "Getting Started"     => "getting-started.md",
        "Types and Functions" => "functions.md",
        "Reference Index"     => "reference.md",
    ],
)

makedocs(;
    format = Documenter.HTML(
        canonical = "http://tkluck.github.io/LandauGinzburgCategories.jl/stable/",
    ),
    config...
)

makedocs(;
    format = LaTeX(
        platform = "docker",
    ),
    config...
)

deploydocs(
    repo   = "github.com/tkluck/LandauGinzburgCategories.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
