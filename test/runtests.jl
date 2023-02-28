using Revise
using SpectralOperators
using Makie
using CairoMakie
using KrylovKit
# using Test




# @testset "SpectralOperators.jl" begin
#     # Write your tests here.
# end


N = 200
Space = (0.,1.)
kind = :Cheb1
SpectDisc = SpectralDiscretization1D(Space, :Cheb1, N)
T = CoeffTransform(SpectDisc)
T⁻¹ = RealTransform(SpectDisc)
∫̂ = IntegrateCoef(SpectDisc)

∫ = f -> 1/(SpectDisc.Space[2] - SpectDisc.Space[1])*(T⁻¹∘∫̂∘T)(f)

function ∫₀ˣ(f)

    ∫f = ∫(f)

    return ∫f .- ∫f[end]

end

x = zeros(N)

for i=1:N

    x[i] = (cos((i-1)*π/(N-1)) +1)/2

end

κ = 1/10
ψₙ = (1 .+ x*κ)
φₙ = exp.(-10*x.^2)
k = 20

K = fₙ -> 2*κ ./ψₙ .* ∫₀ˣ(fₙ)  - k^2* 4*x .* ∫₀ˣ( x./ψₙ.^2 .* ∫₀ˣ(fₙ) )  
b = 4*x .* ∫₀ˣ(φₙ)
A = fₙ -> fₙ + K(fₙ)

v, info = linsolve(A, b)
u = ∫(v)

#plot the solution
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "u")
lines!(ax, x, v, color = :red, label = "u")
lines!(ax, x[1:end-1], v[1:end-1]./x[1:end-1], color = :red, label = "u")
lines!(ax, x, u, color = :blue, label = "φₙ")
fig



function const_matrix()

    M = Array{Float64}(undef, N, N)

    for i=1:N
        eᵢ = zeros(N)
        eᵢ[i] = 1.
        mᵢ = (T∘A)(eᵢ)
        for j=1:N
            M[j,i] = mᵢ[j]
        end

    end

    return M


end

M = const_matrix()

λ, v = eigen(M)



λ_r = real.(λ)
λ_i = imag.(λ)

scatter(λ_r, λ_i)





