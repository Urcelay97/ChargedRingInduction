"""
Module created by Nicolás Urcelay.

For futher information please send a mail: nico_mx4@hotmail.com
"""
module carv

export integTrap, integSimpson38, derivate, newRap, LegendrePol, LegendreDer,
legendreRoots, GaussianWeights, integGaussQaud, GaussTemp, integSubGaussQuad

    #Trapezoidal method
"""
    integTrap(f,a,b,N)

Integrates a function `f` in the interval `[a,b]` by trapezoidal method with `N` subintervals.
"""
    function integTrap(f::Function,a::Real,b::Real,N::Integer)

        dx = (b-a)/(N-1)
        area = (f(a)+f(b))*dx/2
        x = a+dx    

        for _ in 2:N-1
            area+= f(x)*dx
            x += dx
        end    

        return area
    end

    #Simpson's 3/8 method
"""
    integSimpson38(f::Function,a::Real,b::Real,N::Integer)

Integrates a function `f` in the interval `[a,b]` with the Simpson's 3/8 method.

`N` must be odd.

"""
    function integSimpson38(f::Function,a::Real,b::Real,N::Integer)

        dx = (b-a)/(N-1)
        area = (f(a)+f(b))*dx/3
        x = a+dx

        for i in 2:N-1
            
            if i%2 == 0
                area += (4/3)*dx*f(x)
                x += dx 
                continue
            end
            area += (2/3)*dx*f(x)
            x += dx 
        end

        return area
    end

    #Numerical first derivate
"""
    derivate(f::Function, x0::Real, h::Real, mode::String)

Performs de first derivate of `f` evaluated in `x0` with a step size of `h`.

# Parameters
## `f::Function`
Scalar function of one variable.
## `x0::Real`
The point at which the first derivative is found.
## `h::Real`
Step size.
## `modes::String`
Methods to perform the derivate: `right` `left` `center` `extended`.

"""
    function derivate(f::Function,x0::Real, h::Real,mode::String)

        if mode == "right"
            return (f(x0 + h) - f(x0))/h
        
        elseif mode == "left"
            return (f(x0) - f(x0 - h))/h

        elseif mode == "center"
            return (f(x0 + h) - f(x0 - h))/(2*h)

        elseif mode == "extended"
            return (4(f(x0 + h/2) - f(x0 - h/2))/(h) - (f(x0 + h) - f(x0 - h))/(2*h))/3
        end

    end

    #Newton Rapshon
"""
Obtains the first root of a function `f` with the Newton's - Raphson's method given an initial value `x0`.
# Parameters
## `f::Function`
Scalar function of one variable.
## `x0::Real`
Initial value.
## `h::Real`
Step size for the derivate in the method.
## `error::Real`
Precision to achieve.
"""
    function newRap(f::Function,x0::Real,h::Real,error::Real)
        
        x = x0 - f(x0)/derivate(f,x0,h,"extended")

        while abs(x-x0)>=error
            x0 = x
            x = x0 - f(x)/derivate(f,x0,h,"extended")
        end

        return x

    end

    #Legendre Polynomials
"""
    LegendrePol(x::Real,N::Integer)

Evaluates the `N` Lengendre's polynomial in `x`

"""
    function LegendrePol(x::Real,N::Integer)

        sum = 0
        big(N)
        for i in 0:N
            combCoef = factorial(big(N))/(factorial(big(i))*factorial(big(N-i)))
            sum += (combCoef)^2*((x+1)^(N-i))*((x-1)^i)
        end
        
        return sum/2^N
    end

    #First derivate of Lengendre's Polynomial
"""
    LegendreDer(x::Real,N::Integer)

Performs the first derivate of the `N` Legendre's polynomial evaluated in `x`
"""
    function LegendreDer(x::Real,N::Integer)

        if N == 0
            return 0

        elseif N == 1
            return 1

        else
            FirstTerm = N*(x+1)^(N-1)
            LastTerm = N*(x-1)^(N-1)
            sum = FirstTerm + LastTerm

            for i in 1:N-1
                combCoef = factorial(big(N))/(factorial(big(i))*factorial(big(N-i)))
                sum += (combCoef)^2*( (N-i)*(x+1)^(N-i-1)*(x-1)^i + (x+1)^(N-i)*i*(x-1)^(i-1) )
            end
            return sum/2^N
        end
    end

    #Legendre's polynomial roots
"""
    LegendreRoots(N::Integer,error::Real)

Return an array with the `N` roots of the `N` Legendre's Polynomial with a precision of `error`.
"""
    function LegendreRoots(N::Integer,error::Real)
        
        #First iteration
        roots = zeros(N)
        #Newton Raphson with the exact derivate
        x0 = -1
        x = x0 - LegendrePol(x0,N)/LegendreDer(x0,N)
        while abs(x-x0)>=error
            x0 = x
            x = x0 - LegendrePol(x0,N)/LegendreDer(x0,N)
        end
        roots[1] = x

        #Next iteration (Numerical derivate)
        for i in 2:N
            f(t) = LegendrePol(t,N)/prod(t.-roots[1:i-1])
            roots[i] = newRap(f,x0,1e-3,error)
        end

        return roots
    end

    #Gaussian quadrature weights
"""
    GaussianWeights(N::Integer, error::Real)

Calculates the `N` Gaussian Quadrature weights of `N` points with a precision of `error`
"""
    function GaussianWeights(N::Integer, error::Real)
        
        roots = LegendreRoots(N, error)
        weights = zeros(N)

        for (ind, root) in enumerate(roots)
            weights[ind] = 2/( (1-root^2)*(LegendreDer(root,N))^2 )
        end
        return weights
    end
    
    #Gaussian quadrature
"""
    integGaussQaud(f::Function,a::Real,b::Real,N::Integer,error::Real)

Performs the Gaussian Quadrature integration of the function `f` in the interval `[a,b]` with `N` points
given an `error` in the Gaussian weights and Legendre's roots.
"""
    function integGaussQuad(f::Function,a::Real,b::Real,N::Integer,error::Real)

        weights = GaussianWeights(N, error)
        x = LegendreRoots(N, error)
        sum = 0

        for i in 1:N
            sum += weights[i]*f( ((b-a)/2)*x[i] + (a+b)/2 )
        end

        return sum*(b-a)/2
    end

    #Gaussian quadrature without calculation of weights
"""
    GaussTemp(f::Function,a::Real,b::Real,N::Integer,error::Real,w::AbstractArray{Real})

Performs the Gaussian Quadrature integration of the function `f` in the interval `[a,b]` with `N` points
given an `error` in the Legendre's roots.

The values of the Gaussian weights `w` must be provided in an array.
"""
    function GaussTemp(f::Function,a::Real,b::Real,N::Integer,error::Real,w)
        x = LegendreRoots(N, error)
        sum = 0

        for i in 1:N
            sum += w[i]*f( ((b-a)/2)*x[i] + (a+b)/2 )
        end
        return sum*(b-a)/2
    end

    #Gaussian quadrature with subdomains
"""
    integSubGaussQuad(f::Function,a::Real,b::Real,N::Integer,error::Real,n::Integer)

Performs the Gaussian Quadrature integration `n` times in `n-1` subintervals with `N` points
in order to integrate a function `f` in the interval `[a,b]` given an `error` in the Legendre's roots
and Gaussian weights.

"""
    function integSubGaussQuad(f::Function,a::Real,b::Real,N::Integer,error::Real,n::Integer)
        sum = 0
        w = GaussianWeights(N, error)
        dx = (b-a)/(n-1)
        low = a
        for _ in 1:n-1
            sup = low + dx
            sum += GaussTemp(f,low,sup,N,error,w)
            low = sup
        end
        return sum
    end
    
end


                


    




