include("CalcNumCTEC.jl")

f(x,y) = 2x*y + 2x - x^2 - 2y^2 + 72

# f(x,y) = 1 + x + y

a = 0.
b = 8.

c = 0.
d = 6.

I_ap = IntegralDupla_Trapezio(f, (a,b), (c,d))
I_ex = 2544.

println("Integral aproximada: ", I_ap)
println("Integral exata     : ", I_ex)
println("Erro absoluto      : ", abs(I_ap - I_ex))