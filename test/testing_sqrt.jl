using Revise
using SpectralOperators
using Makie
using CairoMakie
using KrylovKit
using LinearAlgebra


N = 500
R = 1.0
Space = (0.,R)
kind = :Cheb1
SpectDisc = SpectralDiscretization1D(Space, :Cheb1, N);
CT = CoeffTransform(SpectDisc);
CT⁻¹ = RealTransform(SpectDisc);
# ∫̂ = IntegrateCoef(SpectDisc)
∫̂l = IntegrateCoefL(SpectDisc);
∫̂r = IntegrateCoefR(SpectDisc);

# ∫ = f -> 1/(SpectDisc.Space[2] - SpectDisc.Space[1])*(T⁻¹∘∫̂∘T)(f)

∫l = f -> 1/(SpectDisc.Space[2] - SpectDisc.Space[1])*(CT⁻¹∘∫̂l∘CT)(f)

∫r = f -> 1/(SpectDisc.Space[2] - SpectDisc.Space[1])*(CT⁻¹∘∫̂r∘CT)(f)


x = zeros(N)


for i=1:N

    x[i] = (cos((i-1)*π/(N-1)) +1)* (SpectDisc.Space[2] -SpectDisc.Space[1] )/2 - SpectDisc.Space[1] 

end


#condition

ns = collect(158:0.25:190)
condition = 0.0*ns

κ = 1/4
ψ = (1 .+ x*κ)
φ = exp.(-x.^2)
# n = 1.
for i ∈ eachindex(ns)
    n = ns[i]
    T = σ -> 1/2*( (x/R).^2 .- 1) .*  ∫l(x.^2 .* σ) + (x.^2)/2 .* ∫r(( (x/R).^2 .- 1).*σ)
    T′ = σ -> x/R^2 .*  ∫l(x.^2 .* σ) + x.* ∫r(( (x/R).^2 .- 1).*σ)

    Pₙ = σ -> σ + 2*(κ ./ ψ ) .* T′(σ) - (4*n^2) * (x ./ ψ) .* T(σ)

    A = zeros(N,N)
    # P̂ₙ = Pₙ∘CT
    for i=1:N
        eᵢ = zeros(N)
        eᵢ[i] = 1
        A[:,i] = (CT⁻¹∘Pₙ∘CT)(eᵢ)
    end
    condition[i] = cond(A)
end

#plot condition number in log scale
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel = "n", ylabel = "condition number", xscale = log10, yscale = log10)
lines!(ax, ns, condition, color = :red, label = "condition number")
fig
#estimate the order of growth of the condition number
(log.(condition)[end] - log.(condition)[end-1])/(log.(ns)[end] - log.(ns)[end-1])



(log.(condition)[end] - log.(condition)[end-1])/(log.(ns)[end] - log.(ns)[end-1])

T = σ -> 1/2*( (x/R).^2 .- 1) .*  ∫l(x.^2 .* σ) + (x.^2)/2 .* ∫r(( (x/R).^2 .- 1).*σ)
T′ = σ -> x/R^2 .*  ∫l(x.^2 .* σ) + x.* ∫r(( (x/R).^2 .- 1).*σ)

Pₙ = σ -> σ + 2*(κ ./ ψ ) .* T′(σ) - (4*n) * (x ./ ψ) .* T(σ)

σ = zeros(N)
σ, info = linsolve(Pₙ , φ)

u = T(σ)
u′ = T′(σ)
#plot results
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "u(x)")
lines!(ax, x, u, color = :red, label = "u(x)")
lines!(ax, x[1:end], u′[1:end], color = :blue, label = "u′(x)")
#axislegend("", position = :lt)
axislegend( position = :lt)
fig


#generate Matrix
#We will generate it in coef Space





   

