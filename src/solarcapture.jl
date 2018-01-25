"""
	solarcapturesimple(sunshine, trajectory, solarpanels)

Calculate the Flux, Power, and Total Energy for all input panels for all input time.
"""
function solarcapturesimple(sunshine, trajectory, solarpanels)

	#initialize
	timestep = sunshine.time[2]-sunshine.time[1]
	panelflux = Array{Float32,2}(length(solarpanels),length(sunshine.time))
	panelpower = Array{Float32,2}(length(solarpanels),length(sunshine.time))
	paneltotalenergy = Array{Float32,2}(length(solarpanels),length(sunshine.time))

	for i=1:length(sunshine.time)

		# --- Get Normalized Solar Vector --- #
		A = (-sunshine.azimuth[sunshine.time[i]] + 90.0)*pi/180.0 #Correct azimuth angle orientation (x-axis = 0 , positive counter-clockwise)
		Z = (sunshine.zenith[sunshine.time[i]])*pi/180.0

		# Convert sun angles to vector
		u = cos(A).*sin(Z)
		v = sin(A).*sin(Z)
		w = cos(Z)

		sunvector = [u,v,w]
		# ---

		for j=1:length(solarpanels)

			# --- Adjust panel normal by aircraft orientation --- #
			# complete trig calculations for brevity
			cphi = cos(-trajectory.phi[i])
			cth = cos(-trajectory.theta[i])
			cpsi = cos(trajectory.psi[i])
			sphi = sin(-trajectory.phi[i])
			sth = sin(-trajectory.theta[i])
			spsi = sin(trajectory.psi[i])

			# 3D Rotation Matrix
			R = [cth*cpsi, -cphi*spsi+sphi*sth*cpsi,  sphi*spsi+cphi*sth*cpsi;
				 cth*spsi,  cphi*cpsi+sphi*sth*spsi, -sphi*cpsi+cphi*sth*spsi;
				-sth, 	    sphi*cth, 				 cphi*cth]

			panelnormal = R*solarpanels.normal[j]

			# --- Get obliquity factor (0-1) --- #
			mu = sum(sunvector.*panelnormal)
			if mu < 0.0
				mu = 0.0
			end

			# --- Calculate Flux, Power, and Total Energy --- #
			panelflux[j,i] = sunshine.flux[sunshine.time[i]].*mu

			panelpower[j,i] = PAR[:etasolar]*solarpanel.area[j]*energy.flux[j,i]

			paneltotalenergy[j,i] = panelpower[j,i]*timestep

		end #for number of panels
	end #for all sunshine.time

	return panelenergy(panelflux,panelpower,paneltotalenergy)

end #solarcapturesimple()

"""
    getsolardata(file::String)
Reads SMARTS data from a file and stores the data as a global sunshine type
named SMARTSdata.
"""
function getsolardata(file::String)
  if isfile(file)
    data = readtable(file, header = true)
  else
    error("SMARTS DATA file not found at: $file")
  end

  # Extract SMARTS Data from file
  timeSMARTSfile = data[:,1] # time in hours
  fluxSMARTSfile = data[:,3] # total solar flux
  azimuthSMARTSfile = data[:,9]*pi/180 # azimuth angle of solar flux
  zenithSMARTSfile = data[:,10]*pi/180 # zenith angle of solar flux

	global SMARTSdata = sunshine(time,azimuth,zenith,flux)
  return SMARTSdata
end
