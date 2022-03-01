module carv

export integTrap, integSimpson38, integGaussQaud, integSubGaussQuad

    #Integral por metodo del trapecio
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

    #Integral mediante Simpson 3/8
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

    #Derivada numérica
    function derivada(f::Function,x0::Real, h::Real,mode::String)

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
    function newRap(f::Function,x0::Real,h::Real,error::Real)
        
        x = x0 - f(x0)/derivada(f,x0,h,"extended")

        while abs(x-x0)>=error
            x0 = x
            x = x0 - f(x)/derivada(f,x0,h,"extended")
        end

        return x

    end

    #Polinomios de Legendre
    function LegendrePol(x::Real,N::Integer)

        sum = 0
        big(N)
        for i in 0:N
            combCoef = factorial(big(N))/(factorial(big(i))*factorial(big(N-i)))
            sum += (combCoef)^2*((x+1)^(N-i))*((x-1)^i)
        end
        
        return sum/2^N
    end

    #Derivada polinomio de Legendre
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

    #Racies de los polinomios de Legendre
    function LegendreRoots(N::Integer,error::Real)
        
        #Realizando la primera iteración
        roots = zeros(N)
        #Newton Raphson pero con la función exacta de la derivada
        x0 = -1
        x = x0 - LegendrePol(x0,N)/LegendreDer(x0,N)
        while abs(x-x0)>=error
            x0 = x
            x = x0 - LegendrePol(x0,N)/LegendreDer(x0,N)
        end
        roots[1] = x

        #Demás iteraciones (Sin la función exacta)
        for i in 2:N
            f(t) = LegendrePol(t,N)/prod(t.-roots[1:i-1])
            roots[i] = newRap(f,x0,1e-3,error)
        end

        return roots
    end

    #Pesos de la cuadratura Gaussiana
    function GaussianWeights(N::Integer, error::Real)
        
        roots = LegendreRoots(N, error)
        weights = zeros(N)

        for (ind, root) in enumerate(roots)
            weights[ind] = 2/( (1-root^2)*(LegendreDer(root,N))^2 )
        end
        return weights
    end
    
    #Cuadratura Gaussiana
    function integGaussQaud(f::Function,a::Real,b::Real,N::Integer,error::Real)

        weights = GaussianWeights(N, error)
        x = LegendreRoots(N, error)
        sum = 0

        for i in 1:N
            sum += weights[i]*f( ((b-a)/2)*x[i] + (a+b)/2 )
        end

        return sum*(b-a)/2
    end

    #Cuadratura de Gauss sin calcular los pesos
    function GaussTemp(f::Function,a::Real,b::Real,N::Integer,error::Real,w )
        x = LegendreRoots(N, error)
        sum = 0

        for i in 1:N
            sum += w[i]*f( ((b-a)/2)*x[i] + (a+b)/2 )
        end
        return sum*(b-a)/2
    end

    #Cuadratura Gaussiana a intervalos

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


                


    




