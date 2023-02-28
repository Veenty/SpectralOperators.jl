module SpectralOperators

# Write your package code here.

export SpectralDiscretization1D, CoeffTransform, RealTransform, IntegrateCoef, InverseTransform



using FFTW
include("Chebyshev.jl")

struct SpectralDiscretization1D

    Space
    kind #:Fourier, :Cheb1, :Cheb2, :Leg
    N  #Number of points
    A #Matrix needed for change of coordinates, can it be FFT, DCT, LegMatrix, etc.
    A⁻¹ #Inverse of A 
  
end


function SpectralDiscretization1D(Space, kind, N)

    if kind == :Cheb1
        DCT = FFTW.plan_r2r(ones(N) ,FFTW.REDFT00, 1, flags = FFTW.MEASURE)
        return SpectralDiscretization1D(Space, kind, N, DCT, DCT)    
        
    elseif kind == :Cheb2
        DCT = FFTW.plan_r2r(ones(N) ,FFTW.REDFT10, 1, flags = FFTW.MEASURE)
        return SpectralDiscretization1D(Space, kind, N, DCT, DCT)
    else

        error("not yet implemented")
    
    end


end





function CoeffTransform(SpectDisc)
  
    if SpectDisc.kind == :Cheb1
        T = fₙ -> 𝘾(fₙ, SpectDisc.A, 1)
    elseif SpectDisc.kind == :Cheb2
        T = fₙ -> 𝘾(fₙ, SpectDisc.A, 2)
    else
        error("Not yet implemented")
    end

    return T

end

function RealTransform(SpectDisc)

    if SpectDisc.kind == :Cheb1
       T⁻¹ = f̂ -> 𝘾⁻¹(f̂, SpectDisc.A⁻¹, 1)

    elseif SpectDisc.kind == :Cheb2
        T⁻¹ = f̂ -> 𝘾⁻¹(f̂, SpectDisc.A⁻¹, 2)
     else
          error("Not yet implemented")
     end

     return T⁻¹

end


function IntegrateCoef(SpectDisc)

    if SpectDisc.kind == :Cheb1 || SpectDisc.kind == :Cheb2
        ∫̂ = f̂ -> IntegrateCoeff(f̂)
    end

    return ∫̂

end


function InverseTransform(f̂, SpectDisc)

    if SpectDisc.kind == :Cheb1
        return 𝘾⁻¹(f̂, SpectDisc.A⁻¹, 1)
    elseif SpectDisc.kind == :Cheb2
        return 𝘾⁻¹(f̂, SpectDisc.A⁻¹, 2)
    else
        error("Edges not yet implemented")
    end

end









end
