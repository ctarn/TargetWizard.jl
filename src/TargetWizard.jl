module TargetWizard

include("TargetSelect.jl")
include("report/BasicAquisitionReport.jl")
include("report/TargetSelectionReport.jl")
include("view/TargetXView.jl")
include("view/TargetXDualView.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()
main_BasicAquisitionReport()::Cint = BasicAquisitionReport.julia_main()
main_TargetSelectionReport()::Cint = TargetSelectionReport.julia_main()
main_TargetXView()::Cint = TargetXView.julia_main()
main_TargetXDualView()::Cint = TargetXDualView.julia_main()

end
