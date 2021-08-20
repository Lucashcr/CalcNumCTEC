include("CalcNumCTEC.jl")

t0 = time()
function Fatorial(n::Int64)
    if n==0 || n==1
        return 1
    else
        return n * Fatorial(n-1)
    end
end
tf = time()
Fatorial(5)
