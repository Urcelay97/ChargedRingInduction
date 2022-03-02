module InductionFunctions

include("NumericalIntegralFunctions.jl")
using .carv

export Bz, By, Bx

#Induction Z component
function Bz(ro,z,phi,a)
    
    function fz(x)
        num = a^2 - a*ro*cos(phi-x)
        den = sqrt(ro^2 + a^2 + z^2 - 2*a*ro*cos(phi-x))^3
        return num/den
    end

    return integSubGaussQuad(fz,0,2pi,4,1e-5,10)
end

#Induction Y component
function By(ro,z,phi,a)

    function fy(x)
        num = z*a*sin(x)
        den = sqrt(ro^2 + a^2 + z^2 - 2*a*ro*cos(phi-x))^3
        return num/den
    end
    return integSubGaussQuad(fy,0,2pi,4,1e-5,10)
end


#Induction X component
function Bx(ro,z,phi,a)

    function fx(x)
        num = z*a*cos(x)
        den = sqrt(ro^2 + a^2 + z^2 - 2*a*ro*cos(phi-x))^3
        return num/den
    end

    return integSubGaussQuad(fx,0,2pi,4,1e-5,10)
end

end
