import PackageCompiler
import TOML

exes = [
    "TargetSelect", "TargetSelectX", "TargetBind",
    "BasicAquisitionReport", "TargetSelectionReport", "CoverageReport", "CoverageReportX",
    "TargetViewX", "TargetDualViewX", "CrossLinkSiteView", "ExhaustiveSearchViewX",
]

cfg = TOML.parsefile("Project.toml")
dir = "tmp/$(Sys.ARCH).$(Sys.iswindows() ? "Windows" : Sys.KERNEL)/$(cfg["name"])"

PackageCompiler.create_app(".", dir; executables=[exe => "main_$(exe)" for exe in exes],
    force=true, include_lazy_artifacts=true, include_transitive_dependencies=true,
)

open(io -> write(io, cfg["version"]), joinpath(dir, "VERSION"); write=true)
