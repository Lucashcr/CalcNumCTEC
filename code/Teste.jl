include("CalcNumCTEC.jl")

T(x,y) = sin(x) + cos(y) + x*y

IntegDupSimpson(T, 0, 8, 0, 6)