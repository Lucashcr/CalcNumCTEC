include("CalcNumCTEC.jl")

# f(x,y) = 2x*y + 2x - x^2 - 2y^2 + 72
# # f(x,y) = 1 + x + y

# a = 0.
# b = 8.

# c = 0.
# d = 6.

# t0 = time()
# I_ap = IntegralDupla_Trapezio(f, (a,b), (c,d))
# t1 = time()
# I_ex = 2816.

# println("Integral aproximada: ", I_ap)
# println("Integral exata.....: ", I_ex)
# println("Erro absoluto......: ", abs(I_ap - I_ex))
# println("Tempo de execução..: ", round(t1-t0, digits=4), "s")

v1 = Vertex2d(1,2)
println(v)