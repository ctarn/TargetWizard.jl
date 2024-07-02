import LibGit2

import Documenter

repo = "github.com/ctarn/TargetWizard.jl.git"

root = "doc"
out = joinpath(root, "tmp", "doc")

rm(out; force=true, recursive=true)
LibGit2.clone("https://$(repo)", out, branch="gh-pages")
rm(joinpath(out, ".git"); force=true, recursive=true)

vs = readdir(joinpath(root, "log")) .|> splitext .|> first .|> VersionNumber
sort!(vs; rev=true)
logs = map(vs) do v in
    "<section>$(read(joinpath(root, "log", "$(v).html"), String))</section>"
end

html = read(joinpath(root, "index.html"), String)
html = replace(html, "{{ release }}" => "<div class=\"release\">$(join(logs))</div>")
open(io -> write(io, html), joinpath(out, "index.html"); write=true)

open(io -> write(io, "targetwizard.ctarn.io"), joinpath(out, "CNAME"); write=true)

for file in ["fig"]
    cp(joinpath(root, file), joinpath(out, file); force=true)
end

Documenter.deploydocs(repo=repo, target=joinpath("..", out), versions=nothing)

pages = [
    "index.md",
    "guide.md",
    "Tutorial" => ["tutorial/prot.md", "tutorial/xl.md"],
    "Manual" => [
        "manual/index.md",
        "Analysis Reports" => [
            "manual/report/BasicAcquisitionReport.md"
        ],
        "Interactive Visualization Views" => [],
    ],
    "dev.md",
]

Documenter.makedocs(; sitename="TargetWizard", build=joinpath("..", out), pages)
Documenter.deploydocs(; repo, target=joinpath("..", out), dirname="doc")
