# ----- Solar Panel Types ----- #
# Define panelgeometry Type
struct panelgeometry
	normal::Array{Float32,1}(3)
	chord::Float32
	span::Float32
	area::Float32
	roll::Float32
	pitch::Float32
	yaw::Float32
	#TODO talk to Nathaniel about what needs to go in here as far as the power stuff is concerned
end

function panelautodef(;
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

# Define panelphysics type
struct panelphysics
	flux::Array{Float32,2}
	power::Array{Float32,2}
	energy::Array{Float32,2}
end #energy type

# ----- Sunshine Types ----- #
# Define sunshine Type
struct sunshine
	time::Array{Float32,1}
	azimuth::Array{Float32,1}
	zenith::Array{Float32,1}
	flux::Array{Float32,1}
end
