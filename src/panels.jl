#Define Panel Type
struct panel
	normal::Array{Float32,1}(3)
	chord::Float32
	span::Float32
	area::Float32
	roll::Float32
	pitch::Float32
	yaw::Float32
	#TODO talk to Nathaniel about what needs to go in here as far as the power stuff is concerned
end

function paneldefault(;
	normal = [0; 0; 1],
	chord = 1.0,
	span = 1.0,
	area = 1.0,
	roll = 0.0,
	pitch = 0.0,
	yaw = 0.0
	)

	return panel(normal, chord, span, area, roll, pitch, yaw)
end
