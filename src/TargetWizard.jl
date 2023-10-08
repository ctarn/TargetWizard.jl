module TargetWizard

include("TargetSelect.jl")
include("report/BasicAquisitionReport.jl")
include("report/TargetSelectionReport.jl")
include("report/XLCoverageReport.jl")
include("view/TargetXView.jl")
include("view/TargetXDualView.jl")
include("view/XSiteView.jl")
include("view/XExhaustiveSearchView.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()
main_BasicAquisitionReport()::Cint = BasicAquisitionReport.julia_main()
main_TargetSelectionReport()::Cint = TargetSelectionReport.julia_main()
main_XLCoverageReport()::Cint = XLCoverageReport.julia_main()
main_TargetXView()::Cint = TargetXView.julia_main()
main_TargetXDualView()::Cint = TargetXDualView.julia_main()
main_XSiteView()::Cint = XSiteView.julia_main()
main_XExhaustiveSearchView()::Cint = XExhaustiveSearchView.julia_main()

end
