import CSV
import GeometricTools
import Solar
using Interpolations

# set time step
timestep = 1/360.0
t = 0.0:timestep:24.0

# --- Get Solar --- #
data = CSV.read("35_-105_18.288_1.609_2016_12_21_-7.csv"; header=1)
solartime = data[:,1]
flux = data[:,3]
azimuth = data[:,9]*pi/180
zenith = data[:,10]*pi/180
solar = Solar.sunshine(solartime, azimuth, zenith, flux)
solar = Solar.interpolatesunshine(t, solar)

# --- Get Trajectory --- #
ntimesteps = length(t)
roll = fill(atan(35.0^2.0/(9.81*3000)), ntimesteps)
pitch = zeros(ntimesteps)
omega = 3000/35.0*1.0/3600.0; #hours per radian
yaw = zeros(ntimesteps); # start facing north
for i = 2:length(yaw)
    yawstep = (t[i]-t[i-1])/omega; #radians per timestep
    yaw[i] = yaw[i-1]+yawstep
end
trajectory = Solar.aircrafttrajectory(t, roll, pitch, yaw)

# --- Get Geometry --- #
npanel = 10
b = 42.0 #span
dihedral = vcat(fill(0.0, npanel), 80.0, 80.0, 80.0)
chord = vcat(2.0, fill(1.4, npanel), 1.0, 0.8, 0.5)
sweep = fill(20.0*pi/180, length(chord)-1)

t = vcat(0.0, 0.03, range(0.03, stop = 1.0, length = npanel)[2:end], 1.03, 1.04, 1.05)
xyz = zeros(Float64,length(t),3)
for i = 2:length(t)
    dt = t[i]-t[i-1]
    dx = dt*sin(sweep[i-1])
    dy = dt*cos(sweep[i-1])
    dz = dt*cos(sweep[i-1])*tan(dihedral[i-1])
    ds = sqrt(dx^2+dy^2+dz^2)
    dx *= dt/ds
    dy *= dt/ds
    dz *= dt/ds
    xyz[i:end,1] .+= dx
    xyz[i:end,2] .+= dy
    xyz[i:end,3] .+= dz
end
xyz *= b/xyz[end,2]

# reflect geometry
chord = vcat(chord[end:-1:1], chord[2:end])
xpos = vcat(xyz[end:-1:1,1], xyz[2:end,1])
ypos = vcat(-xyz[end:-1:1,2], xyz[2:end,2])
zpos = vcat(xyz[end:-1:1,3], xyz[2:end,3])

# leading edge position of each airfoil
airfoil_files = fill("naca2412.dat", length(chord))
Os = [ [xpos[i], ypos[i], zpos[i]] for i in 1:size(airfoil_files)[1]]

# Orientation of chord of each airfoil (yaw, pitch, roll)
orien = fill([0.0, 0.0, 270.0], length(chord))

# parameterization
n_upper = 5 # Number of sections in upper surface of blade
n_lower = 5 # Number of sections in lower surface of blade
r = fill(1.0, length(chord)) # Expansion ratio in both surfaces of each airfoil
sections = fill([(1.0, 10, 1.0, false)], length(chord)-1) # discretization between each airfoil

crosssections = []
for (i,airfoil_file) in enumerate(airfoil_files)
    # Read airfoil file
    out = CSV.read(airfoil_file; delim = " ", header=0)
    x = out[:,1]
    y = out[:,2]

    # Separate upper and lower sides to make the contour injective in x
    leidx = argmin(x)
    upper = [x[leidx:-1:1], y[leidx:-1:1]]
    lower = [x[leidx:end], y[leidx:end]]

    # Parameterize both sides independently
    fun_upper = GeometricTools.parameterize(upper[1], upper[2], zeros(Float64, length(upper[1])); inj_var=1)
    fun_lower = GeometricTools.parameterize(lower[1], lower[2], zeros(Float64, length(lower[1])); inj_var=1)

    # New discretization for both surfaces
    upper_points = GeometricTools.discretize(fun_upper, 0, 1, n_upper, r[1]; central=true)
    lower_points = GeometricTools.discretize(fun_lower, 0, 1, n_lower, r[1]; central=true)

    # Put both surfaces back together from TE over the top and from LE over the bottom.
    reverse!(upper_points)                       # Trailing edge over the top
    new_x = [point[1] for point in upper_points]
    new_y = [point[2] for point in upper_points] # Leading edge over the bottom
    new_x = vcat(new_x, [point[1] for point in lower_points])
    new_y = vcat(new_y, [point[2] for point in lower_points])

    # Scales the airfoil acording to its chord length
    new_x = chord[i]*new_x
    new_y = chord[i]*new_y

    # Reformats into points
    npoints = size(new_x)[1]
    airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]

    # Positions the airfoil along the blade in the right orientation
    Oaxis = GeometricTools.rotation_matrix(orien[i][1], orien[i][2], orien[i][3])
    invOaxis = inv(Oaxis)
    airfoil = GeometricTools.countertransform(airfoil, invOaxis, Os[i])

    push!(crosssections, airfoil)
end

geometry = zeros(length(crosssections),length(crosssections[1]),length(crosssections[1][1]))
for i = 1:length(crosssections)
    for j = 1:length(crosssections[1])
        geometry[i,j,:] = crosssections[i][j][:]
    end
end

# switch to stability frame
geometry[:,:,1] = -geometry[:,:,1]
geometry[:,:,3] = -geometry[:,:,3]

# get solar panels
points, cells, solarpanels = Solar.rectanglemesh(geometry)
vtk_points = [points[i,:] for i = 1:size(points,1)]
vtk_cells = [cells[i].-1 for i = 1:length(cells)]

# get solar capture
shadowed, unshadowed = Solar.solarcapture(solar, trajectory, solarpanels, 0.25)
#
# # output VTK files
# if !ispath("out")
#     mkdir("out")
# end
#
# for i = 1:size(shadowed.energy,1)
#     if all(shadowed.shadow[i,:])
#         continue
#     end
#     data = []
#     push!(data, Dict(
#                 "field_name" => "Shadow",
#                 "field_type" => "scalar",
#                 "field_data" => Float64.(shadowed.shadow[i,:])
#                 )
#     )
#     push!(data, Dict(
#             "field_name" => "Flux",
#             "field_type" => "scalar",
#             "field_data" => shadowed.flux[i,:]
#             )
#     )
#     push!(data, Dict(
#              "field_name" => "Power",
#              "field_type" => "scalar",
#              "field_data" => shadowed.power[i,:]
#              )
#     )
#     push!(data, Dict(
#               "field_name" => "Energy",
#               "field_type" => "scalar",
#               "field_data" => shadowed.energy[i,:]
#               )
#     )
#     push!(data, Dict(
#                "field_name" => "UnshadowedFlux",
#                "field_type" => "scalar",
#                "field_data" => unshadowed.flux[i,:]
#                )
#     )
#     push!(data, Dict(
#                 "field_name" => "UnshadowedPower",
#                 "field_type" => "scalar",
#                 "field_data" => unshadowed.power[i,:]
#                 )
#      )
#     push!(data, Dict(
#              "field_name" => "UnshadowedEnergy",
#              "field_type" => "scalar",
#              "field_data" => unshadowed.energy[i,:]
#              )
#     )
#     # Generates the VTK file
#     GeometricTools.generateVTK("out/solar$i", vtk_points; cells=vtk_cells, cell_data = data)
#  end
