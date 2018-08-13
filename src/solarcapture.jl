"""
    solarcapturesimple(solardata::sunshine, trajectory::aircrafttrajectory,
        solarpanels::Array{panelgeometry,1}, etasolar::Real)

    Calculate the flux, power, and total energy for all input panels for all
    input time.  Assumes solar panels are in the vehicle frame (NED)
"""
function solarcapturesimple(solardata::sunshine, trajectory::aircrafttrajectory,
    solarpanels::Array{panelgeometry,1}, etasolar::Real)

    #initialize
    flux = Array{Float64,2}(length(solardata.time)-1, length(solarpanels))
    power = Array{Float64,2}(length(solardata.time)-1, length(solarpanels))
    energy = Array{Float64,2}(length(solardata.time)-1, length(solarpanels))

    for i=1:length(solardata.time)-1

        # --- Normalized Solar Vector --- #
        A = solardata.azimuth[i]
        Z = solardata.zenith[i]
        u = cos(A).*sin(Z)
        v = sin(A).*sin(Z)
        w = -cos(Z)
        sunvector = vcat(u,v,w)

        # --- Aircraft rotation matrix --- #
        cphi = cos(trajectory.roll[i])
        cth = cos(trajectory.pitch[i])
        cpsi = cos(trajectory.yaw[i])
        sphi = sin(trajectory.roll[i])
        sth = sin(trajectory.pitch[i])
        spsi = sin(trajectory.yaw[i])
        R = [cth*cpsi -cphi*spsi+sphi*sth*cpsi   sphi*spsi+cphi*sth*cpsi;
             cth*spsi  cphi*cpsi+sphi*sth*spsi  -sphi*cpsi+cphi*sth*spsi;
            -sth            sphi*cth                            cphi*cth]

        # --- Find flux on each panel
        for j=1:length(solarpanels)

            # rotate panel normal with aircraft
            panelnormal = R*solarpanels[j].normal

            # get obliquity factor (0-1)
            mu = sum(sunvector.*panelnormal)
            if mu < 0.0
                mu = 0.0
            end

            # flux, power, and total energy
            flux[i,j] = solardata.flux[i].*mu
            power[i,j] = etasolar*solarpanels[j].area*flux[i,j]
            energy[i,j] = power[i,j]*(solardata.time[i+1]-solardata.time[i])
        end
    end
    return panelenergy(flux, power, energy)
end #solarcapturesimple()

"""
    rectanglemesh(crosssections::Array{Float64,3})

    Creates a rectangular mesh.  Returns points, cells, solarpanels.  The first
    two outputs may be used for visualization with VTK.  The final output is used
    for solar capture calculations.
    For crosssections input dimensions are as follows:
        i = crosssections
        j = points
        k = x,y,z
"""
function rectanglemesh(crosssections::Array{Float64,3})

    # get number of panels
    npanels = (size(crosssections,1)-1)*(size(crosssections,2)-1)

    # prepare solar panel output
    solarpanels = Array{panelgeometry,1}(npanels)

    # prepare vtk output
    points = reshape(crosssections, div(length(crosssections),3), 3)
    cells = Array{Array{Int64,1},1}(npanels)

    # use cross product to find panel area and normal vector
    for i = 1:size(crosssections,1)-1 # loop over sections
        for j = 1:size(crosssections,2)-1 # loop over points in each section
            # panel (i,j), (i+1,j), (i+1,j+1), (i,j+1)
            idx = (size(crosssections,2)-1)*(i-1)+j
            # define rectangular cell
            cells[idx] =
                [sub2ind(crosssections, i, j, 1),
                 sub2ind(crosssections, i+1, j, 1),
                 sub2ind(crosssections, i+1, j+1, 1),
                 sub2ind(crosssections, i, j+1, 1)]
            panelnodes = [crosssections[i,j,:], crosssections[i+1,j,:],
                crosssections[i+1,j+1,:], crosssections[i,j+1,:]]
            # use two triangles to get area and normal
            # triangle 1
            ab = crosssections[i+1,j,:]-crosssections[i,j,:]
            ac = crosssections[i+1,j+1,:]-crosssections[i,j,:]
            tmp1 = cross(ab,ac)
            tmp2 = norm(cross(ab,ac))
            panelarea1 = tmp2/2
            panelnormal1 = tmp1/tmp2
            # triangle 2
            ab = ac
            ac = crosssections[i,j+1,:]-crosssections[i,j,:]
            tmp1 = cross(ab,ac)
            tmp2 = norm(cross(ab,ac))
            panelarea2 = tmp2/2
            panelnormal2 = tmp1/tmp2
            # combine to get area and normal
            panelarea = panelarea1+panelarea2
            panelnormal = (panelnormal1+panelnormal2)/
                norm(panelnormal1+panelnormal2)
            # assemble solar panel structure
            solarpanels[idx] = panelgeometry(panelnodes, panelarea, panelnormal)
        end
    end
    return points, cells, solarpanels
end

"""
    solarcapture(solardata::sunshine, trajectory::aircrafttrajectory,
        solarpanels::Array{panelgeometry,1}, etasolar::Real)

Returns (shadowed, unshadowed) `panelenergy` structs
"""
function solarcapture(solardata::sunshine, trajectory::aircrafttrajectory,
    solarpanels::Array{panelgeometry,1}, etasolar::Real)

    # find solar capture with no shadows
    unshadowed = solarcapturesimple(solardata, trajectory, solarpanels, etasolar)

    # no flux means panel is shadowed
    shadow = unshadowed.flux.==0.0

    # check for additional shadowing
    for i = 1:size(shadow, 1)
        getshadowed!(view(shadow,i,:), trajectory.roll[i],
            trajectory.pitch[i], trajectory.yaw[i], solardata.azimuth[i],
            solardata.zenith[i], solarpanels)
    end

    # apply shadow
    shadowed = deepcopy(unshadowed)
    idxshadow = find(shadow)
    shadowed.shadow[idxshadow] .= true
    shadowed.flux[idxshadow] .= 0.0
    shadowed.power[idxshadow] .= 0.0
    shadowed.energy[idxshadow] .= 0.0

    return shadowed, unshadowed
end #solarcapture

"""
    getshadowed!(shadow::Array{Bool,1}, roll::Real, pitch::Real, yaw::Real,
        azimuth::Real, zenith::Real, solarpanels::Array{panelgeometry,1})
    Returns mask of shadowed panels
"""
function getshadowed!(shadow::AbstractArray{Bool,1}, roll::Real, pitch::Real, yaw::Real,
    azimuth::Real, zenith::Real, solarpanels::Array{panelgeometry,1})
    # return if nothing to calculate
    if all(shadow)
        return nothing
    end

    # error checking
    if length(shadow) != length(solarpanels)
        error("Length of `shadow` must match length of `solarpanels`")
    end

    # get number of panels
    npanels = length(solarpanels)

    # get roll, pitch, yaw rotation
    cphi = cos(roll)
    cth = cos(pitch)
    cpsi = cos(yaw)
    sphi = sin(roll)
    sth = sin(pitch)
    spsi = sin(yaw)
    R = [cth*cpsi -cphi*spsi+sphi*sth*cpsi   sphi*spsi+cphi*sth*cpsi;
         cth*spsi  cphi*cpsi+sphi*sth*spsi  -sphi*cpsi+cphi*sth*spsi;
        -sth            sphi*cth                            cphi*cth]

    # get transformation from vehicle to solar frame
    A = azimuth
    Z = zenith
    sA = sin(A)
    cA = cos(A)
    sZ = sin(Z)
    cZ = cos(Z)
    T  =[sZ*cA sZ*sA -cZ; -sA cA 0; cZ*cA cZ*sA sZ]

    # combine rotation and transformation
    TR = T*R

    # apply rotation and transformation to panels
    sundistance = Array{Float64,1}(npanels)
    panelcentroids2D = Array{Array{Float64,1},1}(npanels)
    panelnodes2D = Array{Array{Array{Float64,1},1},1}(npanels)
    @inbounds for i in eachindex(solarpanels)
        # transform panel coordinates
        panelnodes2D[i] = Array{Array{Float64,1},1}(length(solarpanels[i].nodes))
        for j in eachindex(solarpanels[i].nodes)
            # extract 2D information, discard distance to sun
            panelnodes2D[i][j] = TR[2:3,1:3]*solarpanels[i].nodes[j]
        end
        # transform panel centroids
        panelcentroid = TR*solarpanels[i].centroid
        # split into relative distance to sun and 2D node location
        sundistance[i] = -panelcentroid[1]
        panelcentroids2D[i] = panelcentroid[2:3]
    end

    # scale entire 2D system to 1.0-2.0 for use with GeometricalPredicates package
    maxval = 0.0
    @inbounds for i = 1:length(panelnodes2D)
        for j = 1:length(panelnodes2D[i])
            newmaxval = max(maximum(panelnodes2D[i][j]),abs(minimum(panelnodes2D[i][j])))
            if newmaxval > maxval
                maxval = newmaxval
            end
        end
    end
    factor = 1/(2*maxval)
    offset = 1.5

    @inbounds for i in eachindex(panelnodes2D)
        for j in eachindex(panelnodes2D[i])
            panelnodes2D[i][j] .*= factor
            panelnodes2D[i][j] .+= offset
        end
        panelcentroids2D[i] .*= factor
        panelcentroids2D[i] .+= offset
    end

    # assemble 2D panels as viewed from the sun
    panels2D = Array{GeometricalPredicates.Polygon2D{GeometricalPredicates.Point2D},1}(npanels)
    @inbounds for i in eachindex(panelnodes2D)
        points2D = Array{GeometricalPredicates.Point2D,1}(length(panelnodes2D[i]))
        for j in eachindex(points2D)
            points2D[j] = Point(panelnodes2D[i][j][1], panelnodes2D[i][j][2])
        end
        panels2D[i] = Polygon(points2D...)
    end

    # find shadowed panels
    @inbounds for i = 1:npanels # shadowing panels
        if shadow[i] == false # only unshadowed panels may participate
            for j = 1:npanels # shadowed panels
                if (shadow[j] == false) && (sundistance[j] < sundistance[i]) && (i!=j)
                    # then check if the panel's centroid is shadowed
                    shadow[i] = inpolygon(panels2D[j],
                        Point(panelcentroids2D[i][1], panelcentroids2D[i][2]))
                end
            end
        end
    end
    return nothing
end

"""
    interpolatesunshine(time::AbstractArray{<:Real,1}, solardata::sunshine)
    Interpolates solardata struct data to fit a new array of time values.
"""
function interpolatesunshine(time::AbstractArray{<:Real,1}, solardata::sunshine)
    azimuthitp = interpolate((solardata.time,),solardata.azimuth,Gridded(Linear()))
    zenithitp = interpolate((solardata.time,),solardata.zenith,Gridded(Linear()))
    fluxitp = interpolate((solardata.time,),solardata.flux,Gridded(Linear()))
    return sunshine(time,azimuthitp[time],zenithitp[time],fluxitp[time])
end
