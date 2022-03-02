include("InductionFieldComponents.jl")
using .InductionFunctions
using Plots

#Constants
(XMIN,XMAX) = (-1,1)
(YMIN,YMAX) = (-1,1)
N_POINTS_X = 15
N_POINTS_Y = 15
STEP_LENGTH_X = (XMAX-XMIN)/(N_POINTS_X-1)
STEP_LENGTH_Y = (YMAX-YMIN)/(N_POINTS_Y-1)
RADIUS = 1
Z = -1 #Value of Z which plane XY is ploted
MU = 1 #Magnetic permeability
I = 1 #Current

#X and Y range
x_interval = XMIN:STEP_LENGTH_X:XMAX
y_interval = YMIN:STEP_LENGTH_Y:YMAX

#X,Y,ρ and ϕ meshgrids
coordinates_x = [x for x in x_interval, y in y_interval]' #Transpose in order to create a X meshgrid
coordinates_y = [x for x in x_interval, y in y_interval]
coordinate_ro = sqrt.(coordinates_x.^2 + coordinates_y.^2)
coordinate_phi = atan.(coordinates_y,coordinates_x)

#Calculating X and Y components of induction
induction_x = MU*I/4pi .* Bx.(coordinate_ro,Z,coordinate_phi,RADIUS)
induction_y = MU*I/4pi .* By.(coordinate_ro,Z,coordinate_phi,RADIUS)

#Reshaping meshgrids into vectors to put them into quiver function
x = reshape(coordinates_x,(1,length(coordinates_x)))[1:length(coordinates_x)]
y = reshape(coordinates_y,(1,length(coordinates_y)))[1:length(coordinates_y)]

#Creating a LinRange of Z values to plot B(0,0,z)
z = LinRange(-10,10,1000)

#Reshaping X,Y and Z inductions
induction_x = reshape(induction_x,(1,length(induction_x)))[1:length(induction_x)]
induction_y = reshape(induction_y,(1,length(induction_y)))[1:length(induction_y)]
induction_z = MU*I/4pi .* Bz.(0,z,0,RADIUS)

#Normalizing the lengths of the vectors
#r = @. sqrt(induction_x^2+induction_y^2)

#Plotting
a = quiver(x,y,quiver=(induction_x,induction_y),aspect_ratio=:equal,c=:viridis)
b = plot(z,induction_z,xlabel="a")
plot(a,b,xlabel=["X" "Z"],ylabel=["Y" "Bz"],title=["B(x,y,$Z)" "B(0,0,$Z)"],label=["" ""],dpi=300)