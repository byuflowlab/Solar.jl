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
