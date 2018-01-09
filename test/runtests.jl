using Solar
using Base.Test

tests = ["paneltests", "windtests"]

for t in tests
	include("$(t).jl")
end
