function trial_data()
  flux=zeros(300)
  ψ=zeros(300)
  ϕ=zeros(300)
  α=zeros(300)
  azimuth=zeros(300)
  zenith=zeros(300)
  time=zeros(300)

  flux[1:300]=1376
  zenith[1:300]=20

  for i=1:300
    azimuth[i]=30+i
    time[i]=i
  end

  data=Dict(
  :flux=>flux,
  :ψ=>ψ, #Yaw Angle
  :ϕ=>ϕ, #Bank Angle
  :α=>α,  #Roll angle
  :zenith=>zenith,
  :azimuth=>azimuth,
  :time=>time
  )

  return data
end
