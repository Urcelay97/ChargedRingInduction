#MODULOS

#Modulo con las funciones para la integración numérica
include("integral.jl")
using .carv
#Libreria para graficar
using PlotlyJS

#Main program

#Parametros
a = 1
mu = 1
I = 1

#Definiendo las funciones a integrar

#Bz
function Bz(ro,z,phi,a)
    
    function fz(x)
        num = a^2 - a*ro*cos(phi-x)
        den = sqrt(ro^2 + a^2 + z^2 - 2*a*ro*cos(phi-x))^3
        return num/den
    end

    return integSubGaussQuad(fz,0,2pi,4,1e-5,10)
end

#By
function By(ro,z,phi,a)

    function fy(x)
        num = z*a*sin(x)
        den = sqrt(ro^2 + a^2 + z^2 - 2*a*ro*cos(phi-x))^3
        return num/den
    end

    return integSubGaussQuad(fy,0,2pi,4,1e-5,10)
end

#Bx
function Bx(ro,z,phi,a)

    function fx(x)
        num = z*a*cos(x)
        den = sqrt(ro^2 + a^2 + z^2 - 2*a*ro*cos(phi-x))^3
        return num/den
    end

    return integSubGaussQuad(fx,0,2pi,4,1e-5,10)
end

#Definiendo las posiciones a graficar

#Definiendo los limites
x = LinRange(-1.5,1.5,15)
y = LinRange(-1.5,1.5,15)
z = LinRange(-1.5,1.5,15)

#Creando el plano xy
X = zeros(length(x),length(y))
Y = zeros(length(x),length(y))

for i in 1:length(x)
    Y[i,:] = y
end

for i in 1:length(y)
    X[:,i] = x
end

XX = zeros(length(x),length(y),length(z))
YY = zeros(length(x),length(y),length(z))
ZZ = zeros(length(x),length(y),length(z))
#Para X y Y
for i in 1:length(z)
    XX[:,:,i] = X
    YY[:,:,i] = Y
end
#Para Z
for i in 1:length(x)
    for j in 1:length(y)
        ZZ[i,j,:] = z
    end
end
#Convirtiendo a arrays de 1 dimensión
XX = reshape(XX,(1,length(XX)))
YY = reshape(YY,(1,length(YY)))
ZZ = reshape(ZZ,(1,length(ZZ)))
#Creando un array para ro y phi
ro = @. sqrt(XX^2 + YY^2)
phi = @. atan(YY/XX)
#Calculando las inducciones
bbx = @. Bx(ro,ZZ,phi,a)
bby = @. By(ro,ZZ,phi,a)
bbz = @. Bz(ro,ZZ,phi,a)
#Para que se vean iguales
r = @. sqrt(bbx^2 + bby^2 + bbz^2)
bbbx = @. bbx/r
bbby = @. bby/r
bbbz = @. bbz/r
#Graficando

circX = [a*cos(t) for t in 0:2pi/100:2pi]
circY = [a*sin(t) for t in 0:2pi/100:2pi]

circulo = scatter3d(
    x=circX, y=circY,z=zeros(100),mode="lines",marker=attr(color="#1f77b4", size=20, symbol="circle",
    line=attr(color="rgb(0,0,0)", width=5)),
line=attr(color="#1f77b4", width=1))

vectores = cone(
    x=XX[1,:],
    y=YY[1,:],
    z=ZZ[1,:],
    u=bbbx[1,:],
    v=bbby[1,:],
    w=bbbz[1,:],
    sizemode="equal",
    sizeref=1,
    anchor="tip"
    )

layout = Layout(autosize=false, width=500, height=500,
    margin=attr(l=0, r=0, b=0, t=65))

plot([circulo,vectores],layout)


Bx(1,1,1,1)

