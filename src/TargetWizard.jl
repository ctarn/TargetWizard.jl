module TargetWizard

include("TargetSelect.jl")
include("AquisitionReport.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()
main_AquisitionReport()::Cint = AquisitionReport.julia_main()

end
