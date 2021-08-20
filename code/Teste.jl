include("CalcNumCTEC.jl")

# T(x,y) = sin(x) + cos(y) + x*y

# IntegDupSimpson(T, 0, 8, 0, 6)
T(x,y,z) = sin(x) + cos(y) + z + sqrt(x*y*z)

t0 = time()
n = 1000
I = IntegralTripla_Simpson(T,0,1,0,1,0,1,nx=n,ny=n,nz=n)
tf = time()
Δt = round(tf-t0)
println("Valor da integral: $I\nTempo de Execução: $Δt s")