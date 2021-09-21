################################
# BIBLIOTECA JULIA CalcNumCTEC #
################################

using LinearAlgebra
using Printf
using Base

function Derivada_MDFCentrada(f::Function, x₀::Number, i::Int; h::Number=10^-4, n::Int=9)
    x = x₀-((n-1)/2)*h:h:x₀+((n-1)/2)*h
    
    A = zeros(n,n)
    for i in range(1,stop=n)
        for j in range(1,stop=n)
            A[i,j] = x[i]^(n-j)
        end
    end

    coefs = A \ f.(x)

    function dif_f_interp(x_in)
        res = 0
        for k in 1:n-i
            res += coefs[k] * ((factorial(n-k))/(factorial(n-k-i))) * x_in^(n-k-i)
        end
        return res
    end

    return dif_f_interp(x₀)
end

# Cálculo de zeros de funções ------------------------------------------------------------------------------

function Zeros_Bissecao(f::Function, a::Number, b::Number; 
    tol::Number=10^-12, klim::Number=10^6)
    if f(a) * f(b) < 0
        k = 0
        x = (a + b) / 2
        
        while abs(f(x)) > tol
            if k == klim
                @assert false "ERRO: O método não convergiu para $klim iterações..."
            else
                k += 1
            end
        
            x = (a + b) / 2
            if f(x) * f(a) < 0 
                b = x
            else
                a = x 
            end
        end

        return x
    else
        @assert false "\n\nERRO NA FUNÇÃO ZeroBissecao:\nOs valores da função nos extremos do intervalo dado\nprecisam ter valores de sinais opostos entre si\npara que a convergência do método seja garantida..."
    end
end

function Zeros_Cordas(f::Function, a::Number, b::Number; 
    tol::Float64=10^-12, klim::Int64=10^6)
    if f(a) * f(b) < 0
        k = 0
        x = (a * f(b) - b * f(a)) / (f(b) - f(a))
        
        while abs(f(x)) > tol
            if k == klim
                @assert false "ERRO: O método não convergiu para $klim iterações..."
            else
                k += 1
            end
        
            x = (a + b) / 2
            if f(x) * f(a) < 0
                b = x 
            else
                a = x 
            end
        end

        return x
    else
        @assert false "ERRO NA FUNÇÃO ZeroCordas:\nOs valores da função nos extremos do intervalo dado",
        "\nprecisam ter valores de sinais opostos entre si",
        "\npara que a convergência do método seja garantida..."
    end
end

function Zeros_NR(f::Function, x::Number; 
    tol::Float64=10^-12, klim::Int64=10^6)
    k = 0
    while abs(f(x)) > tol
        if k == klim
            @assert false "ERRO NA FUNÇÃO ZeroNR:\nO método não convergiu para $klim iterações..."
        else
            k += 1
        end
        
        if abs(Derivada_MDFCentrada(f, x, 1, n=5)) < tol
            @assert false "ERRO NA FUNÇÃO ZeroNR:\nO ponto x=$x possui derivada nula..."
        end
        x = x - (f(x) / Derivada_MDFCentrada(f, x, 1, n=5))
    end

    return x
end

# Resolução de sistemas de equações lineares --------------------------------------------------------------

function SEL_EliminGauss(A::Matrix, b::Vector; tol::Float64=10^-12)
    (m,n) = size(A)
    
    if m==n && m == length(b)
        if abs(det(A)) < tol
            @assert false "ERRO na função EliminGauss:\nA matriz é singular..."
        else
            function test(M)
                aux = 0
                for j in range(1, stop = n)
                    for i in range(1, stop = n)
                        if i > j 
                            aux += abs(M[i,j])
                        end
                    end
                end
                return aux
            end

            AG = copy(A)
            bG = copy(b)
            
            for j in 1:n
                i = j
                while A[j,j] == 0
                    i+=1

                    aux_A = AG[j,:]
                    AG[j,:] = AG[i,:]
                    AG[i,:] = aux_A

                    aux_b = b[j]
                    bG[j] = bG[i]
                    bG[i] = aux_b
                    if (i,j) == (m,n)
                        @assert false "ERRO na função EliminGauss:\nA matriz é singular..."
                    end
                end
            end
        
            i = 1
            j = 1
            k = 0
            while test(AG) > tol
                k += 1
                if i == n
                    j += 1
                    i = j
                end
                i += 1
        
                bG[i] -= bG[j] * (AG[i,j] / AG[j,j])
                AG[i,:] -= AG[j,:] * (AG[i,j] / AG[j,j])
            end
        
            x = zeros(n)
            for k in n:-1:1
                x[k] = (bG[k] - dot(AG[k,:], x)) / AG[k,k]
            end
        
            return x
        end
    else
        @assert false "ERRO NA FUNÇÃO ElimGauss:\nInconsistência nas dimensões da matriz e do vetor..."
    end
end

function SEL_CritConverg(A::Matrix, modo::Char)    
    (m,n) = size(A)

    if m == n
        if uppercase(modo) == 'L'
            beta=zeros(n)
            for i in range(1, stop = n)
                for j in range(1, stop = n)
                    beta[i] += abs(A[i,j]) 
                end
                beta[i] /= abs(A[i,i])
            end
            return max(beta...)
        elseif uppercase(modo) == 'C'
            beta=zeros(n)
            for j in range(1, stop = n)
                for i in range(1, stop = n)
                    beta[j] += abs(A[i,j]) 
                end
                beta[j] /= abs(A[j,j])
            end
            return max(beta...)
        elseif uppercase(modo) == 'S'
            beta = ones(n)
            for i in range(1, stop = n)
                for j in range(1, stop = n)
                    beta[i] += abs(A[i,j]) * beta[j]
                end
                beta[i] -= 1
                beta[i] /= abs(A[i,i])
            end
            return max(beta...)
        end
    else
        @assert false "\n\nERRO NA FUNÇÃO CritConverg:\nA matriz dada não é quadrada..."
    end
end

function SEL_GaussJacobi(A::Matrix, b::Vector; x_in::Vector=zeros(length(b)), 
    tol::Float64=1e-12, klim::Int=10^6)
    # Rseolve o sistema linear A*x=b através de iterações pelo Método da Gauss-Jacobi.
    # Tolerância padrão: 10^(-12)
    # Métrica de erro: eₐ = ||x_out - x_in|| / ||x_out||
    # Retorno da (x_out, err)

    (m,n) = size(A)

    if m==n && m == length(b)
        if abs(det(A)) < tol
            @assert false "ERRO na função GaussJacobi:\nA matriz é singular..."
        else
            for j in 1:n
                i = j
                while A[j,j] == 0
                    i+=1

                    aux_A = A[j,:]
                    A[j,:] = A[i,:]
                    A[i,:] = aux_A

                    aux_b = b[j]
                    b[j] = b[i]
                    b[i] = aux_b
                    if (i,j) == (m,n)
                        @assert false "ERRO na função EliminGauss:\nA matriz é singular..."
                    end
                end
            end

            function err(x1, x2)
                return norm(x2 - x1) / norm(x2)
            end
            
            k = 1
            x_out = zeros(n)
            for i in range(1, stop = n)
                x_out[i] = (b[i] - dot(A[i,1:i-1], x_in[1:i-1]) - dot(A[i,i+1:n], x_in[i+1:n])) / A[i,i]
            end
            
            while err(x_in, x_out) > tol
                if k == klim
                    @assert false "ERRO NA FUNÇÃO GaussJacobi:\nO método não convergiu para $klim iterações\ne tolerância de $tol..."
                else
                    k += 1
                end
            
                x_in = copy(x_out)
            
                for i in range(1, stop = n)
                    x_out[i] = (b[i] - dot(A[i,1:i-1], x_in[1:i-1]) - dot(A[i,i+1:n], x_in[i+1:n])) / A[i,i]
                end
            end
            
            return x_out
        end
    else
        @assert false "ERRO NA FUNÇÃO GaussJacobi:\nInconsistência nas dimensões da matriz e do vetor..."
    end
end

function SEL_GaussSeidel(A::Matrix, b::Vector; x_in::Vector=zeros(length(b)), 
    tol::Float64=1e-12, klim::Int=10^6)
    # Rseolve o sistema linear A*x=b através de iterações pelo Método da Gauss-Seidel.
    # Tolerância padrão: 10^(-12)
    # Métrica de erro: eₐ = ||x_out - x_in|| / ||x_out||
    # Retorno da (x_out, err)

    (m,n) = size(A)

    if m==n && m == length(b)
        if abs(det(A)) < tol
            @assert false "\n\nERRO na função GaussJacobi:\nA matriz é singular...\n\n"
        else
            for j in 1:n
                i = j
                while A[j,j] == 0
                    i+=1

                    aux_A = A[j,:]
                    A[j,:] = A[i,:]
                    A[i,:] = aux_A

                    aux_b = b[j]
                    b[j] = b[i]
                    b[i] = aux_b
                    if (i,j) == (m,n)
                        @assert false "ERRO na função GaussSeidel:\nA matriz é singular..."
                    end
                end
            end

            function err(x1, x2)
                return norm(x2 - x1) / norm(x2)
            end
            
            k = 1
            x_out = zeros(n)
            for i in range(1, stop = n)
                x_out[i] = (b[i] - dot(A[i,1:i-1], x_out[1:i-1]) - dot(A[i,i+1:end], x_in[i+1:end])) / A[i,i]
            end
            
            while err(x_in, x_out) > tol
                if k == klim
                    @assert false "ERRO NA FUNÇÃO GaussSeidel:\nO método não convergiu para $klim iterações\ne tolerância de $tol..."
                else
                    k += 1
                end
            
                x_in = copy(x_out)
            
                for i in range(1, stop = n)
                    x_out[i] = (b[i] - dot(A[i,1:i-1], x_out[1:i-1]) - dot(A[i,i+1:end], x_in[i+1:end])) / A[i,i]
                end
            end
            
            return x_out
        end
    else
        printstyled("ERRO NA FUNÇÃO GaussJacobi:\nInconsistência nas dimensões da matriz e do vetor...\n\n", color=:bold)
    end
end

# Resolução de sistemas de equações não lineares ---------------------------------------------------------

function SENL_NRSecante(F::Function, x_in::Vector;
    tol::Float64=1e-6, klim::Int64=10^3)
    # Rseolve o sistema não linear F(x)=0 através de iterações pelo Método de Newton-Raphson Secante.
    # Ou seja, a Jacobiana é calculada de maneira aproximada em função da tolerância.
    # Tolerância padrão: 10^(-12)

    function J_aprox(x, δ=tol)
        n = length(x)
        Jap = zeros(n,n)
    
        for j in range(1, stop = n)
            δx = zeros(n)
            δx[j] = δ * x[j]
            for i in range(1, stop = n)
                Jap[i,j] = (F(x + δx)[i] - F(x)[i]) / δx[j]
            end
        end
    
        return Jap
    end

    if length(F(x_in)) == length(x_in)
        err(x_in, x_out) = norm(x_out - x_in) / norm(x_out)

        k = 1
        x_out = x_in - J_aprox(x_in) \ F(x_in)
        
        while err(x_in, x_out) > tol
            if k==klim
                @assert false "\n\nERRO na função NR_Secante:\nO método não convergiu para $klim iterações\ne tolerância de $tol...\n\n"
            else
                k+=1
            end
        
            x_in = copy(x_out)
            x_out = x_in - J_aprox(x_in) \ F(x_in)
        end
        
        return x_out
    else
        printstyled("\n\nERRO na função NR_Secante:\nOs vetores possuem dimensões diferentes...\n\n", color=:bold)
        return nothing
    end
end

# Interpolação e ajuste ------------------------------------------------------------------------------

function Interpolacao_Polinomial(x::Vector, y::Vector)
    # Calcula e retorna o polinômio interpolador para os dados de entrada via polinômios de Lagrange.
    # O grau do polinômio é determinado pelas dimensões dos vetores x e y.

    n = Int(length(x))

    if n==length(y)
        function f_interpoladora(x_in)
            res = 0
            for i in range(1, stop=n)
                aux = 1
                for j in range(1, stop=n)
                    if i != j
                        aux *= (x_in - x[j])/(x[i] - x[j])
                    end
                end
                res += y[i] * aux
            end 
            return res
        end 

        return f_interpoladora
    else
        @assert false "ERRO NA FUNÇÃO Interpolacao_Polinomial:\nOs vetores de entrada possuem dimensões diferentes..."
    end
end

function Interpolacao_Polinomial(f::Function, x::Vector)
    # Calcula e retorna o polinômio interpolador para os dados de entrada via polinômios de Lagrange.
    # O grau do polinômio é determinado pela dimensão do vetor x.

    n = Int(length(x))

    function f_interpoladora(x_in)
        res = 0
        for i in range(1, stop=n)
            aux = 1
            for j in range(1, stop=n)
                if i != j
                    aux *= (x_in - x[j])/(x[i] - x[j])
                end
            end
            res += f(x[i]) * aux
        end 
        return res
    end
    
    return f_interpoladora
end

function Interpolacao_Bilinear(x::Vector, y::Vector, z::Vector)
    if length(x) == length(y) == length(z) == 4
        M = [ones(n) x y x.*y]
        b = M \ z
        g(x,y) = dot(b,[1 x y x*y])
        return g
    else
        @assert false "\n\nERRO NA FUNÇÃO Interpolacao_Bilinear:\nOs vetores de entrada precisam ter dimensões de 4 elementos..."
    end
end

function Ajuste_MinimosQuadrados(x::Vector, y::Vector; grau::Int64=1)
    # Calcula e retorna o polinômio ajustado para os dados de entrada via Método dos Mínimos Quadrados.
    # O grau padrão do polinômio é 1.
    # Retorno da função: (polinomio, r²)
    if length(x) != length(y)
        @assert false "ERRO NA FUNÇÃO MinimosQuadrados:\nOs vetores de entrada possuem dimensões diferentes..."
    else
        A = zeros(grau+1, grau+1)
        for i in 0:grau
            for j in 0:grau
                A[i+1,j+1] = sum(x.^(i+j))
            end
        end
    
        b = zeros(grau+1)
        for i in 0:grau
            b[i+1] = sum(y.*(x.^i))
        end 
        
        coefs = A \ b
    end

    function f_ap(x_in)
        xᵏ = []
        for i in 0:grau
            push!(xᵏ, x_in^i)
        end
        return dot(coefs, xᵏ)
    end

    Sr = sum((y - f_ap.(x)).^2)
    ym = sum(y)/length(y)
    St = sum((y .- ym).^2)
    r² = 1 - (Sr/St)
    
    return f_ap, r²
end

# Integração Numérica ------------------------------------------------------------------

# OBS: COMPLETAR O ERRO COM A DIFERENCIAÇÃO NUMÉRICA

function Integral_Trapezio(f::Function, x::Tuple; tol::Float64 = 1e-12)
    if length(x) == 2
        a = x[1]
        b = x[2]
        c = (a+b)/2

        T(a,b) = (b-a) * (f(a) + f(b)) / 2

        if T(a,b) - T(a,c) - T(c,b) < tol
            return T(a,b)
        else
            return Integral_Trapezio(f, (a,c), tol=tol) + Integral_Trapezio(f, (c,b),tol=tol)
        end
    else
        @assert false "ERRO NA FUNÇÃO Intergral_Trapezio:\nO intervalo precisa ter apenas 2 valores..."
    end
end

function Integral_Trapezio(x::Vector, y::Vector)
    n = length(x)

    if n == length(y)
        I = 0
        for i in 1:n-1
            I += ((x[i+1]-x[i])/2)*(y[i] + y[i+1])
        end

        #err = abs(-(b-a)^4 * (IntegSimpson(dif_f4, a, b, n=10^3))[1] / (180*n^4))

        return I#, err
    else
        @assert false "ERRO NA FUNÇÃO Integral_Trapezio:\nOs vetores possuem dimensões diferentes..."
    end
end

function IntegralDupla_Trapezio(f::Function, x::Tuple, y::Tuple; tol::Float64=1e-12)
    if length(x) == length(y) == 2
        x₀, x₁ = x
        y₀, y₁ = y
        xₘ = (x₀ + x₁) / 2
        yₘ = (y₀ + y₁) / 2

        function T(x,y)
            return Integral_Trapezio([y...], [Integral_Trapezio([x...], f.([x...], [y[1],y[1]])),Integral_Trapezio([x...], f.([x...], [y[2],y[2]]))])
        end

        if T(x,y) - T((x₀,xₘ),y) - T((xₘ,x₁),y) < tol
            if T(x,y) - T(x,(y₀,yₘ)) - T(x,(yₘ,y₁)) < tol
                println("teste1")
                return T(x,y)
            else
                println("teste2")
                return IntegralDupla_Trapezio(f, x, (y₀,yₘ), tol=tol) + IntegralDupla_Trapezio(f, x, (yₘ,y₁), tol=tol)
            end
        elseif T(x,y) - T(x,(y₀,yₘ)) - T(x,(yₘ,y₁)) < tol
            println("teste3")
            return IntegralDupla_Trapezio(f, (x₀,xₘ), y, tol=tol) + IntegralDupla_Trapezio(f, (xₘ,x₁), y, tol=tol)
        else
            println("teste4")
            return IntegralDupla_Trapezio(f, (x₀,xₘ), (y₀,yₘ), tol=tol) + IntegralDupla_Trapezio(f, (xₘ,x₁), (y₀,yₘ), tol=tol) + IntegralDupla_Trapezio(f, (x₀,xₘ), y, tol=tol) + IntegralDupla_Trapezio(f, (xₘ,x₁), y, tol=tol)
        end
    else
        @assert false "ERRO NA FUNÇÃO IntegralDupla_Trapezio:\nOs intervalos precisam ter apenas 2 valores..."
    end
end

function Integral_Simpson(f::Function, x::Tuple; tol::Float64 = 1e-12)
    if length(x) == 2
        a = x[1]
        b = x[2]
        c = (a+b)/2
    
        S(a,b) = (b-a) * (f(a) + 4*f((a+b)/2) + f(b)) / 6

        if S(a,b) - S(a,c) - S(c,b) < tol
            return S(a,b)
        else
            return Integral_Simpson(f, (a,c), tol=tol) + Integral_Simpson(f, (c,b), tol=tol)
        end
    else
        @assert false "ERRO NA FUNÇÃO Integral_Simpson:\nO intervalo precisa ter apenas 2 valores..."
    end
end

function Integral_Simpson(x::Vector, y::Vector)
    n = length(x)

    if n == length(y)
        if n%2==1
            I = 0
            for i in 1:2:n-2
                I += ((x[i+2]-x[i])/6)*(y[i] + 4*y[i+1] + y[i+2])
            end

            return I
        else
            @assert false "ERRO NA FUNÇÃO Integral_Simpson:\nO número de elementos dos vetores deve ser ímpar..."
        end
    else
        @assert false "ERRO NA FUNÇÃO IntegSimpson:\nOs vetores possuem dimensões diferentes..."
    end
end

# function IntegralDupla_Simpson(f::Function, x1::Number, x2::Number, 
#     y1::Number, y2::Number; nx::Number=(b-a)*10^3, ny::Number=(d-c)*10^3)

#     hx = (x2-x1)/2nx
#     hy = (y2-y1)/2ny

#     M = zeros(2nx+1, 2ny+1)
#     for (j,y) in enumerate(y1:hy:y2)
#         for (i,x) in enumerate(x1:hx:x2)
#             M[i,j] = f(x,y)
#         end
#     end

#     V = zeros(2ny+1)
#     for k in 1:2ny+1
#         V[k] += Integral_Simpson(Vector(x1:hx:x2), M[:,k]) 
#     end

#     return Integral_Simpson(Vector(y1:hy:y2), V) 
# end

# function IntegralDupla_Simpson(f::Function, x₀::Number, x₁::Number, y₀::Number, y₁::Number; tol::Float64=1e-12)
#     S(a,b) = (b-a) * (f(a) + 4*f((a+b)/2) + f(b)) / 6

#     xₘ = (x₀ + x₁) / 2
#     yₘ = (y₀ + y₁) / 2

#     if 

#     else
        
#     end
# end

function IntegralTripla_Simpson(f::Function, x1::Number, x2::Number, y1::Number, y2::Number, 
    z1::Number, z2::Number; nx::Number=(x2-x1)*100, ny::Number=(y2-y1)*100, nz::Number=(z2-z1)*100)

    hx = (x2-x1)/2nx
    hy = (y2-y1)/2ny
    hz = (z2-z1)/2nz    

    I = []

    M = Matrix{Number}(undef, 2ny+1, 2nz+1)
    for (k,z) in enumerate(z1:hz:z2)
        for (j,y) in enumerate(y1:hy:y2)
            for (i,x) in enumerate(x1:hx:x2)
                M[i,j] = f(x,y,z)
            end
        end
        # println(M)
        V = zeros(2nx+1)
        for k in 1:2nx+1
            V[k] = Integral_Simpson(Vector(x1:hx:x2), M[:,k]) 
        end
        # println(V)
        push!(I, Integral_Simpson(Vector(y1:hy:y2), V))
    end
    res = Integral_Simpson(Vector(z1:hz:z2), I)
    return res

    # M = Matrix{Number}(undef, 2nx+1, 2ny+1)
    # for (j,y) in enumerate(y1:hy:y2)
    #     for (i,x) in enumerate(x1:hx:x2)
    #         M[i,j] = Integral_Simpson(Vector(x1:hx:x2), S[i,j,:]) 
    #     end
    # end

    # V = zeros(2ny+1)
    # for k in 1:2ny+1
    #     V[k] = Integral_Simpson(Vector(x1:hx:x2), M[:,k]) 
    # end

    # return Integral_Simpson(Vector(y1:hy:y2), V)
end

function IntegralTripla_Trapezio(f::Function, x1::Number, x2::Number, y1::Number, y2::Number, 
    z1::Number, z2::Number; nx::Int=10^3, ny::Int=10^3, nz::Int=10^3)

    hx = (x2-x1)/nx
    hy = (y2-y1)/ny
    hz = (z2-z1)/nz    

    I = []

    M = Matrix{Number}(undef, ny+1, nz+1)
    for (k,z) in enumerate(z1:hz:z2)
        for (j,y) in enumerate(y1:hy:y2)
            for (i,x) in enumerate(x1:hx:x2)
                M[i,j] = f(x,y,z)
            end
        end
        # println(M)
        V = zeros(nx+1)
        for k in 1:nx+1
            V[k] = Integral_Trapezio(Vector(x1:hx:x2), M[:,k]) 
        end
        # println(V)
        push!(I, Integral_Trapezio(Vector(y1:hy:y2), V))
    end
    return Integral_Trapezio(Vector(z1:hz:z2), I)

    # M = Matrix{Number}(undef, 2nx+1, 2ny+1)
    # for (j,y) in enumerate(y1:hy:y2)
    #     for (i,x) in enumerate(x1:hx:x2)
    #         M[i,j] = Integral_Simpson(Vector(x1:hx:x2), S[i,j,:]) 
    #     end
    # end

    # V = zeros(2ny+1)
    # for k in 1:2ny+1
    #     V[k] = Integral_Simpson(Vector(x1:hx:x2), M[:,k]) 
    # end

    # return Integral_Simpson(Vector(y1:hy:y2), V)
end

# Geometria -----------------------------------------------------------------------------------

struct Vertex2d
    x::Number
    y::Number

    function Vertex2d(x::Number, y::Number)
        
    end

    function Base.print(v::Vertex2d)
        print("(", v.x, ", ", v.y, ")")
    end
    function Base.println(v::Vertex2d)
        println("(", v.x, ", ", v.y, ")")
    end
end



struct Vertex3d
    x::Number
    y::Number
    z::Number

    function Vertex3d(v::Vertex2d,z::Number)
        return Vertex3d(v.x, v.y, z)
    end

    function Base.print(v::Vertex3d)
        print("(", v.x, ", ", v.y, ", ", v.z, ")")
    end
    function Base.println(v::Vertex3d)
        println("(", v.x, ", ", v.y, ", ", v.z, ")")
    end
end

"CalcNumCTEC incluída com sucesso!"

