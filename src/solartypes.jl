# ----- Trajectory Types ------ #
"""
   aircrafttrajectory

Holds aircraft orientation information for each time step
"""
struct aircrafttrajectory
	roll::Array{Float64,1}
	pitch::Array{Float64,1}
	yaw::Array{Float64,1}
end

# ----- Solar Panel Types ----- #
"
	panelgeometry

Describe the geometry/orientation of a solar panel array.
"
struct panelgeometry
	normal::Array{Float64,1}
	chord::Float64
	span::Float64
	area::Float64
	roll::Float64
	pitch::Float64
	yaw::Float64
end

"""
	panelgeometry(;chord = 1.0, span = 1.0, roll = 0.0, pitch = 0.0, yaw = 0.0))

Define the panelgeometry immutable type; calculate and fill in the area and normal fields from the required inputs for the user.

See documentation for the panelgeometry immutable type for information on inputs.
"""
function panelgeometry( #TODO find dimensions of single solar cell and put those as defaults
	chord = 1.0,
	span  = 1.0,
	roll  = 0.0,
	pitch = 0.0,
	yaw   = 0.0)

	area = chord*span

	# Calculate panel normal in inertial frame
	c1 = cos(-roll)
	c2 = cos(-pitch)
	c3 = cos(yaw)
	s1 = sin(-roll)
	s2 = sin(-pitch)
	s3 = sin(yaw)

	n1 = c1.*s2.*s3 - c3.*s1
	n2 = c1.*c3.*s2 + s1.*s3
	n3 = c1.*c2

	normal = [n1, n2, n3]

	return panelgeometry(normal, chord, span, area, roll, pitch, yaw)

end

"""
	panelparameters

"""
struct panelparameters
	#TODO talk to Nathaniel about what needs to go in here as far as the power stuff is concerned
end #panelphysicsin

"""
	panelenergy

Define flux, power, and total energy arrays for a solar panel network.
"""
struct panelenergy
	flux::Array{Float64,2}
	power::Array{Float64,2}
	totalenergy::Array{Float64,2}
end #panelphysicsout

# ----- Sunshine Types ----- #
"""
	sunshine

Define time, sun angles (azimuth, zenith), and solar flux arrays.
"""
struct sunshine
	time::Array{Float64,1}
	azimuth::Array{Float64,1}
	zenith::Array{Float64,1}
	flux::Array{Float64,1}
end
