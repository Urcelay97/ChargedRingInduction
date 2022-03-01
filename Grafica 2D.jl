include("integral.jl")
using .carv
using Plots
using ProgressMeter
using FFTW

(xmin,xmax) = (-1,1)
(ymin,ymax) = (-1,1)
N_POINTS_X = 15
N_POINTS_Y = 15
TIME_STEP_LENGTH_X = (xmax-xmin)/(N_POINTS_X-1)
TIME_STEP_LENGTH_Y = (ymax-ymin)/(N_POINTS_Y-1)
RADIUS = 1
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

x_interval = xmin:TIME_STEP_LENGTH_X:xmax
y_interval = ymin:TIME_STEP_LENGTH_Y:ymax

coordinates_x = [x for x in x_interval, y in y_interval]'
coordinates_y = [x for x in x_interval, y in y_interval]
coordinate_ro = sqrt.(coordinates_x.^2 + coordinates_y.^2)
coordinate_phi = atan.(coordinates_y,coordinates_x)

induction_x = Bx.(coordinate_ro,-1,coordinate_phi,RADIUS)
induction_y = By.(coordinate_ro,-1,coordinate_phi,RADIUS)


x = reshape(coordinates_x,(1,length(coordinates_x)))[1:length(coordinates_x)]
y = reshape(coordinates_y,(1,length(coordinates_y)))[1:length(coordinates_y)]
z = LinRange(-10,10,1000)


induction_x = reshape(induction_x,(1,length(induction_x)))[1:length(induction_x)]
induction_y = reshape(induction_y,(1,length(induction_y)))[1:length(induction_y)]
induction_z = Bz.(0,z,0,RADIUS)
r = @. sqrt(induction_x^2+induction_y^2)
quiver(x,y,quiver=(0.01.*induction_x./r,0.01.*induction_y./r),aspect_ratio=:equal)
plot(z,induction_z,dpi=300)