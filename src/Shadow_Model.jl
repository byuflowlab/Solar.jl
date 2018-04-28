using ArrayFire
function Shadow_Model(Mesh,Conditions)
    Precision=Mesh[:precision]
    for i=1:length(Conditions[:flux])
        #Load all input angles
        zenith=-Conditions[:zenith][i]
        azimuth=Conditions[:azimuth][i]
        ψ=Conditions[:ψ][i]*180/pi
        α=-Conditions[:α][i]*180/pi
        ϕ=Conditions[:ϕ][i]*180/pi

        #Calculate rotation matrices about each asxis
        Rz=[cosd(ψ) -sind(ψ) 0;
          sind(ψ) cosd(ψ)  0;
          0        0       1]
        Ry=[cosd(α)  0 sind(α);
          0        1       0;
          -sind(α) 0 cosd(α)]
        Rx=[1 0       0       ;
          0 cosd(ϕ) -sind(ϕ);
          0 sind(ϕ) cosd(ϕ) ]

        Rzsun=[cosd(azimuth) -sind(azimuth) 0;
               sind(azimuth)  cosd(azimuth) 0;
               0              0             1]
        Rysun=[cosd(zenith)   0             sind(zenith);
               0              1             0           ;
              -sind(zenith)   0             cosd(zenith)]

        Data=Array{Float32}(9,length(Mesh[:x]))

        Rotated_Coordinates=transpose(round.(([Mesh[:x][:] Mesh[:y][:] Mesh[:z][:]]*Rz*Ry*Rx*Rzsun*Rysun)/Precision))
        Rotated_Normals=transpose([Mesh[:normal][:i][:] Mesh[:normal][:j][:] Mesh[:normal][:k][:]]*Rzsun*Rysun*Rx*Ry*Rz)

        Data[1,:]=Rotated_Coordinates[1,:]
        Data[2,:]=Rotated_Coordinates[2,:]
        Data[3,:]=Rotated_Coordinates[3,:]
        Data[4,:]=Rotated_Normals[1,:]
        Data[5,:]=Rotated_Normals[2,:]
        Data[6,:]=Rotated_Normals[3,:]
        Data[7,:]=Mesh[:area]
        Data[8,:]=Mesh[:panel]

        Atup = reinterpret(NTuple{3, Float64}, Rotated_Coordinates, (size(Rotated_Coordinates, 2),));
        algorithim=sortperm(Atup; alg=QuickSort, by=x->(x[2],x[3],x[1]),rev=true)

        Data[9,algorithim[1]]=true
        current_y=Data[2,algorithim[1]]
        current_z=Data[3,algorithim[1]]
        top_x=Data[1,algorithim[1]]
        for h=2:length(Data[2,:])
            if Data[2,algorithim[h]]!=current_y||Data[3,algorithim[h]]!=current_z||Data[1,algorithim[h]]==top_x
                current_y=Data[2,algorithim[h]]
                current_z=Data[3,algorithim[h]]
                top_x=Data[1,algorithim[h]]
                Data[9,algorithim[h]]=true
            else
                Data[9,algorithim[h]]=false
            end
        end

        paneled_area=0
        paneled_area_lit=0
        flux=0
        for h=1:length(Data[1,:])
            if Data[8,algorithim[h]]==true
                paneled_area+=Data[7,algorithim[h]]
                if Data[9,algorithim[h]]==true
                    paneled_area_lit+=Data[7,algorithim[h]]
                    percent_perpendicular=dot([Data[4,algorithim[h]],Data[5,algorithim[h]],Data[6,algorithim[h]]],[1,0,0])
                    if percent_perpendicular<0
                        #print("Error")
                    else
                        flux+=percent_perpendicular*Conditions[:flux][i]
                    end
                end
            end
        end
        print("Paneled Area=$paneled_area m^2\nPaneled Area Lit=$paneled_area_lit m^2\nFlux=$flux W\n")

        Dummy=Dict()
        Dummy[:x]=Data[1,:]*Precision
        Dummy[:y]=Data[2,:]*Precision
        Dummy[:z]=Data[3,:]*Precision
        Dummy[:panel]=Data[8,:]
        Dummy[:lit]=Data[9,:]
        return (Dummy)
    end
end
