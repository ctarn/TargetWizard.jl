import LibGit2

import Documenter

import TargetWizard

repo = "github.com/ctarn/TargetWizard.jl.git"

root = "doc"
tmp = "tmp"
out = joinpath(root, tmp)
dst = "doc"

@info "cleaning $(out)"
rm(out; force=true, recursive=true)
mkpath(out)
try
    LibGit2.clone("https://$(repo)", out, branch="gh-pages")
catch e
    @warn e
end
rm(joinpath(out, ".git"); force=true, recursive=true)

vs = readdir(joinpath(root, "log")) .|> splitext .|> first .|> VersionNumber
sort!(vs; rev=true)
logs = map(vs) do v in
    "<section>$(read(joinpath(root, "log", "$(v).html"), String))</section>"
end
html = read(joinpath(root, "index.html"), String)
html = replace(html, "{{ release }}" => "<div class=\"release\">$(join(logs))</div>")
open(io -> write(io, html), joinpath(out, "index.html"); write=true)

for file in ["fig"]
    cp(joinpath(root, file), joinpath(out, file); force=true)
end

Documenter.deploydocs(; repo, target=tmp, versions=nothing, cname="targetwizard.ctarn.io")

pages = [
    "index.md",
    "guide.md",
    "Tutorial" => ["tutorial/prot.md", "tutorial/xl.md"],
    "Manual" => ["manual/index.md", "manual/report.md", "manual/view.md"],
    "dev.md",
]

Documenter.makedocs(; sitename="TargetWizard", build=joinpath(tmp, dst), pages)
Documenter.deploydocs(; repo, target=joinpath(tmp, dst), dirname=dst)
