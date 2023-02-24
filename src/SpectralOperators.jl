module SpectralOperators

# Write your package code here.


using FFTW

struct SpectralDiscretization1D

    Space
    kind #:Fourier, :Cheb1, :Cheb2, :Leg
    N  #Number of points
    A #Matrix needed for change of coordinates, can it be FFT, DCT, LegMatrix, etc.
    A⁻¹ #Inverse of A 
    
    if Type == :Cheb1
        DCT = FFTW.plan_r2r(ones(N) ,FFTW.REDFT00, 1, flags = FFTW.MEASURE)
        SpectralDiscretization1D(Space, kind, N) = new(Space, kind, N, DCT, DCT)    
        A = FFTW.plan_r2r(ones(N) ,FFTW.REDFT00, 1, flags = FFTW.MEASURE)
        A⁻¹ = A

    elseif Type == :Cheb2
        DCT = FFTW.plan_r2r(ones(N) ,FFTW.REDFT10, 1, flags = FFTW.MEASURE)
        SpectralDiscretization1D(Space, kind, N) = new(Space, kind, N, DCT, DCT)


    else

        error("Edges not yet implemented")
    end




end


include("Chebyshev.jl")


function Transform(fₙ, SpectDisc)

    if SpectDisc.kind == :Cheb1
        return 𝘾(fₙ, SpectDisc.A, 1)
    elseif SpectDisc.kind == :Cheb2
        return 𝘾(fₙ, SpectDisc.A, 2)
    else
        error("Edges not yet implemented")
    end

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
