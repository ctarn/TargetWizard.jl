module TargetWizard

include("TargetSelect.jl")
include("TargetSelectX.jl")
include("TargetBind.jl")
include("report/BasicAquisitionReport.jl")
include("report/TargetSelectionReport.jl")
include("report/CoverageReport.jl")
include("report/CoverageReportX.jl")
include("view/TargetViewX.jl")
include("view/TargetDualViewX.jl")
include("view/CrossLinkSiteView.jl")
include("view/ExhaustiveSearchViewX.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()
main_TargetSelectX()::Cint = TargetSelectX.julia_main()
main_TargetBind()::Cint = TargetBind.julia_main()
main_BasicAquisitionReport()::Cint = BasicAquisitionReport.julia_main()
main_TargetSelectionReport()::Cint = TargetSelectionReport.julia_main()
main_CoverageReport()::Cint = CoverageReport.julia_main()
main_CoverageReportX()::Cint = CoverageReportX.julia_main()
main_TargetViewX()::Cint = TargetViewX.julia_main()
main_TargetDualViewX()::Cint = TargetDualViewX.julia_main()
main_CrossLinkSiteView()::Cint = CrossLinkSiteView.julia_main()
main_ExhaustiveSearchViewX()::Cint = ExhaustiveSearchViewX.julia_main()

end
