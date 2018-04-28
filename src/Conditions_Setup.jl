function create_tellemetry()
    ϕ=0#Roll angle radians
    α=0#Pitch angle in radians
    ψ=0#Yaw angle in radians
    time=0
    azimuth=102#In degrees
    zenith=10#In degrees
    flux=1
    return (aircrafttrajectory(ϕ,α,ψ),sunshine(time,azimuth,zenith,flux))
end
