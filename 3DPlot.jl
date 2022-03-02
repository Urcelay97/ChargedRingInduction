include("InductionFieldComponents.jl")
using .InductionFunctions
using PlotlyJS
using ProgressMeter

#CONSTANTS
(XMIN, XMAX) = (-1, 1)
(YMIN, YMAX) = (-1, 1)
(ZMIN, ZMAX) = (-1, 1)
N_POINTS_X = 10
N_POINTS_Y = 10
N_POINTS_Z = 10
STEP_LENGTH_X = (XMAX - XMIN) / (N_POINTS_X - 1)
STEP_LENGTH_Y = (YMAX - YMIN) / (N_POINTS_Y - 1)
STEP_LENGTH_Z = (ZMAX - ZMIN) / (N_POINTS_Z - 1)
RADIUS = 1
MU = 1 #Magnetic permeability
I = 1 #Current

# X,Y and Z range
x_interval = XMIN:STEP_LENGTH_X:XMAX
y_interval = YMIN:STEP_LENGTH_Y:YMAX
z_interval = ZMIN:STEP_LENGTH_Z:ZMAX

# Creating containers for the variables
coordinates_x = zeros(N_POINTS_X*N_POINTS_Y*N_POINTS_Z)
coordinates_y = zeros(N_POINTS_X*N_POINTS_Y*N_POINTS_Z)
coordinates_z = zeros(N_POINTS_X*N_POINTS_Y*N_POINTS_Z)
# Induction
induction_x = zeros(N_POINTS_X*N_POINTS_Y*N_POINTS_Z)
induction_y = zeros(N_POINTS_X*N_POINTS_Y*N_POINTS_Z)
induction_z = zeros(N_POINTS_X*N_POINTS_Y*N_POINTS_Z)

# Obtaining the Field components
@showprogress for (k,z) in enumerate(z_interval)
    for (j,y) in enumerate(y_interval)
        for (i,x) in enumerate(x_interval)
            #Coordinates
            coordinates_x[i + (j-1)*N_POINTS_X + (k-1)*N_POINTS_X*N_POINTS_Y] = x
            coordinates_y[i + (j-1)*N_POINTS_X + (k-1)*N_POINTS_X*N_POINTS_Y] = y
            coordinates_z[i + (j-1)*N_POINTS_X + (k-1)*N_POINTS_X*N_POINTS_Y] = z
            #Induction
            induction_x[i + (j-1)*N_POINTS_X + (k-1)*N_POINTS_X*N_POINTS_Y] = MU*I/4pi*Bx(sqrt(x^2 + y^2),z,atan(y,x),RADIUS)
            induction_y[i + (j-1)*N_POINTS_X + (k-1)*N_POINTS_X*N_POINTS_Y] = MU*I/4pi*By(sqrt(x^2 + y^2),z,atan(y,x),RADIUS)
            induction_z[i + (j-1)*N_POINTS_X + (k-1)*N_POINTS_X*N_POINTS_Y] = MU*I/4pi*Bz(sqrt(x^2 + y^2),z,atan(y,x),RADIUS)
        end
    end
end

#In order to normalize the vectors (Optional)

"""
r = @. sqrt(induction_x^2 + induction_y^2 + induction_z^2)
induction_x = @. induction_x / r
induction_y = @. induction_y / r
induction_z = @. induction_z / r
r = r./maximum(r)
"""

#Plotting
ringX = [RADIUS * cos(t) for t in 0:2pi/100:2pi]
ringY = [RADIUS * sin(t) for t in 0:2pi/100:2pi]

#Ring in the plane (x,y,0)
ring = scatter3d(x = ringX, y = ringY, z = zeros(100), mode = "lines")

vectors = cone(
    x = coordinates_x,
    y = coordinates_y,
    z = coordinates_z,
    u = induction_x,
    v = induction_y,
    w = induction_z,
    sizeref = 1,
    sizemode = "scaled",
    anchor = "tip"
)

layout = Layout(autosize = false, width = 500, height = 500,aspectratio=(x=1, y=1, z=0.8), camera_eye=attr(x=-1.57, y=1.36, z=0.58))

plot([ring, vectors], layout)



