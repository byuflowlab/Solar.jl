"""
    solarcapturesimple(solardata::sunshine, trajectory::aircrafttrajectory,
	  solarpanels::panelgeometry,etasolar::Float64)

Calculate the flux, power, and total energy for all input panels for all input time.
"""
function solarcapturesimple(solardata::sunshine,trajectory::aircrafttrajectory,
	solarpanels::Array{panelgeometry,1},etasolar::Float64)

	#initialize
	panelflux = Array{Float64,2}(length(solarpanels),length(solardata.time))
	panelpower = Array{Float64,2}(length(solarpanels),length(solardata.time))
	paneltotalenergy = Array{Float64,2}(length(solarpanels),length(solardata.time))

	for i=1:length(solardata.time)-1

		# --- Get Normalized Solar Vector --- #
		A = -solardata.azimuth[i] + 90.0*pi/180.0 #Correct azimuth angle orientation (x-axis = 0 , positive counter-clockwise)
		Z = solardata.zenith[i]

		# Convert sun angles to vector
		u = cos(A).*sin(Z)
		v = sin(A).*sin(Z)
		w = cos(Z)

		sunvector = [u,v,w]
		# ---

		for j=1:length(solarpanels)

			# --- Adjust panel normal by aircraft orientation --- #
			# complete trig calculations for brevity
			cphi = cos(-trajectory.roll[i])
			cth = cos(-trajectory.pitch[i])
			cpsi = cos(trajectory.yaw[i])
			sphi = sin(-trajectory.roll[i])
			sth = sin(-trajectory.pitch[i])
			spsi = sin(trajectory.yaw[i])

			# 3D Rotation Matrix
			R = [cth*cpsi -cphi*spsi+sphi*sth*cpsi   sphi*spsi+cphi*sth*cpsi;
				   cth*spsi  cphi*cpsi+sphi*sth*spsi  -sphi*cpsi+cphi*sth*spsi;
				  -sth 	     sphi*cth 				         cphi*cth               ]

			panelnormal = R*solarpanels[j].normal

			# --- Get obliquity factor (0-1) --- #
			mu = sum(sunvector.*panelnormal)
			if mu < 0.0
				mu = 0.0
			end

			# --- Calculate Flux, Power, and Total Energy --- #
			panelflux[j,i] = solardata.flux[i].*mu

			panelpower[j,i] = etasolar*solarpanels[j].area*panelflux[j,i]

			paneltotalenergy[j,i] = panelpower[j,i]*(solardata.time[i+1]-solardata.time[i])

		end #for number of panels
	end #for all solardata.time

	return panelenergy(panelflux,panelpower,paneltotalenergy)

end #solarcapturesimple()

"""
    getsunshine(file::String)
Reads SMARTS data from a file and returns a solardata struct
"""
function getsunshine(file::String)
  if isfile(file)
    data = CSV.read(file, header = 1)
  else
    error("SMARTS DATA file not found at: $file")
  end

  # Extract SMARTS Data from file
  time = data[:,1] # time in hours
  flux = data[:,3] # total solar flux
  azimuth = data[:,9]*pi/180 # azimuth angle of solar flux
  zenith = data[:,10]*pi/180 # zenith angle of solar flux

	return sunshine(time,azimuth,zenith,flux)
end

"""
    interpolatesunshine(solardata::sunshine,time::Array{Float64,1})
Interpolates solardata struct data to fit a new array of time values.
"""
function interpolatesunshine(solardata::sunshine,time::Array{Float64,1})
  azimuthitp = interpolate((solardata.time,),solardata.azimuth,Gridded(Linear()))
  zenithitp = interpolate((solardata.time,),solardata.zenith,Gridded(Linear()))
	fluxitp = interpolate((solardata.time,),solardata.flux,Gridded(Linear()))
	return sunshine(time,azimuthitp[time],zenithitp[time],fluxitp[time])
end
