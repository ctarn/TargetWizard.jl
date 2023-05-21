module TargetWizard

include("TargetSelect.jl")
include("report/BasicAquisitionReport.jl")
include("report/TargetSelectionReport.jl")
include("TargetView.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()
main_BasicAquisitionReport()::Cint = BasicAquisitionReport.julia_main()
main_TargetSelectionReport()::Cint = TargetSelectionReport.julia_main()
main_TargetView()::Cint = TargetView.julia_main()

end
