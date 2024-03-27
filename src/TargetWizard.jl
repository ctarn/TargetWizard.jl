module TargetWizard

include("TargetSelect.jl")
include("TargetSelectX.jl")
include("TargetBind.jl")
include("report/BasicAquisitionReport.jl")
include("report/TargetSelectionReport.jl")
include("report/CoverageReport.jl")
include("report/CoverageReportX.jl")
include("view/TargetXView.jl")
include("view/TargetXDualView.jl")
include("view/CrossLinkSiteView.jl")
include("view/XExhaustiveSearchView.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()
main_TargetSelectX()::Cint = TargetSelectX.julia_main()
main_TargetBind()::Cint = TargetBind.julia_main()
main_BasicAquisitionReport()::Cint = BasicAquisitionReport.julia_main()
main_TargetSelectionReport()::Cint = TargetSelectionReport.julia_main()
main_CoverageReport()::Cint = CoverageReport.julia_main()
main_CoverageReportX()::Cint = CoverageReportX.julia_main()
main_TargetXView()::Cint = TargetXView.julia_main()
main_TargetXDualView()::Cint = TargetXDualView.julia_main()
main_CrossLinkSiteView()::Cint = CrossLinkSiteView.julia_main()
main_XExhaustiveSearchView()::Cint = XExhaustiveSearchView.julia_main()

end
