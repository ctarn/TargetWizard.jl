module TargetWizard

include("TargetSelect.jl")

main_TargetSelect()::Cint = TargetSelect.julia_main()

end
