import Documenter

import TargetWizard

repo = "github.com/ctarn/TargetWizard.jl.git"
root = "doc"
tmp = "tmp"
dst = "doc"

out = joinpath(root, tmp) |> mkpath
html = read(joinpath(root, "index.html"), String)
logs = map((VersionNumber∘first∘splitext).(readdir(joinpath(root, "log"))) |> sort) do v
    "<section>$(read(joinpath(root, "log", "$(v).html"), String))</section>"
end |> reverse |> join
html = replace(html, "{{ release }}" => "<div class=\"release\">$(logs)</div>")
open(io -> write(io, html), joinpath(out, "index.html"); write=true)
for file in ["fig"]
    cp(joinpath(root, file), joinpath(out, file); force=true)
end
Documenter.deploydocs(; repo, target=tmp, versions=nothing, cname="targetwizard.ctarn.io")

pages = [
    "index.md",
    "guide.md",
    "Tutorial" => joinpath.("tutorial", ["prot.md", "xl.md"]),
    "Manual" => joinpath.("manual", ["index.md", "report.md", "view.md"]),
    "dev.md",
]
Documenter.makedocs(; sitename="TargetWizard", build=joinpath(tmp, dst), pages)
Documenter.deploydocs(; repo, target=joinpath(tmp, dst), dirname=dst)
