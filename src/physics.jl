"""
	energysimple(sunshine, trajectory, solarpanels)
"""
function solarcapturesimple(sunshine, trajectory, solarpanels)

	#initialize
	timestep = sunshine.time[2]-sunshine.time[1]
	panelflux = Array{Float32,2}(length(solarpanels),length(sunshine.time))
	panelpower = Array{Float32,2}(length(solarpanels),length(sunshine.time))
	panelenergy = Array{Float32,2}(length(solarpanels),length(sunshine.time))

	for i=1:length(sunshine.time)
		# --- Get Normalized Solar Vector --- #
		A = (-sunshine.azimuth[sunshine.time[i]] + 90.0)*pi/180.0 #Correct azimuth angle orientation (x-axis = 0 , positive counter-clockwise)
		Z = (sunshine.zenith[sunshine.time[i]])*pi/180.0

		# Convert sun angles to vector
		u = cos(A).*sin(Z)
		v = sin(A).*sin(Z)
		w = cos(Z)

		sun_vector = [u,v,w]

		for j=1:length(solarpanels)
			# --- Adjust panel normal by aircraft orientation --- #
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

			#TODO take care of this outside to definition of panels. don't assume symmetric in this function.
			# if k==2
			# 	mirrorplanenormal = [0.0,1.0,0.0] #mirror about the xz plane
			# 	panelnormal = panelnormal - 2*(dot(panelnormal,mirrorplanenormal))*mirrorplanenormal
			# end #if other side

			# --- Get obliquity factor (0-1) --- #
			mu = sum(sun_vector.*panelnormal)
			# for sections facing away from sun, set mu to zero
			if mu < 0.0
				mu = 0.0
			end

			# --- Adjust Interpolated Flux by Obliquity Factor --- #
			panelflux[j,i] = sunshine.flux[sunshine.time[i]].*mu

			# --- Calculate panel Power --- #
			panelpower[j,i] = PAR[:etasolar]*solarpanel.area[j]*energy.flux[j,i]

			# --- Calculate panel Energy --- #
			panelenergy[j,i] = panelpower[j,i]*timestep

		end #for number of panels
	end #for all sunshine.time

	return panelphysics(panelflux,panelpower,panelenergy)

end #solarcapturesimple()
