function surface_generator(aircraft,max_seperation=.01) #A function that accepts aircraft geometry as an input and outputs points to approximate the plane surfaces

  #Initialze the output
  meshx=Array{Float64,1}(0)
  meshy=Array{Float64,1}(0)
  meshz=Array{Float64,1}(0)
  meshi=Array{Float64,1}(0)
  meshj=Array{Float64,1}(0)
  meshk=Array{Float64,1}(0)
  mesharea=Array{Float64,1}(0)
  meshpaneled=Array{Bool}(0)

  #Initialize a few variables for calculating distance between points
  distance=0
  first_point=Array{Float32}(3)

  for i=1:length(aircraft.components)#Reads through every component of the plane
    if aircraft.components[i].is=="wing"#Protocol for wing type components

      #Generates the first foil of the wing based on whether it's a vertical or horizontal wing
      if aircraft.components[i].vertical==false
        start_foil=foil_generator(aircraft.components[i].foil[1],aircraft.components[i].chord[1],0,aircraft.components[i].twist[1],aircraft.components[i].LEx[1],aircraft.components[i].LEy[1],aircraft.components[i].LEz[1])
      else
        start_foil=foil_generator(aircraft.components[i].foil[1],aircraft.components[i].chord[1],pi/2,aircraft.components[i].twist[1],aircraft.components[i].LEx[1],aircraft.components[i].LEy[1],aircraft.components[i].LEz[1])
      end

      for j=2:length(aircraft.components[i].chord)#Reads through each wing cross section
        first_foil=start_foil#Sets the first foil of the subsection
        section_width=sqrt((aircraft.components[i].LEx[j]-aircraft.components[i].LEx[j-1])^2+(aircraft.components[i].LEy[j]-aircraft.components[i].LEy[j-1])^2+(aircraft.components[i].LEz[j]-aircraft.components[i].LEz[j-1])^2)
        n=round(section_width/max_seperation+.5)
        end_foil=foil_generator(aircraft.components[i].foil[j],aircraft.components[i].chord[j],aircraft.components[i].dihedral[j-1],aircraft.components[i].twist[j],aircraft.components[i].LEx[j],aircraft.components[i].LEy[j],aircraft.components[i].LEz[j])
        for k=1:n#Cycles through a number of subsections equal to n
          #Intitializes coordinate arrays
          x=Array{Float32}(0)
          y=Array{Float32}(0)
          z=Array{Float32}(0)

          for l=1:length(end_foil[1])#Interpolates end of subsection using the start and end of the section
            push!(x,start_foil[1][l]+(end_foil[1][l]-start_foil[1][l])*k/n)
            push!(y,start_foil[2][l]+(end_foil[2][l]-start_foil[2][l])*k/n)
            push!(z,start_foil[3][l]+(end_foil[3][l]-start_foil[3][l])*k/n)
          end
          second_foil=[x,y,z]

          for l=1:length(end_foil[1])-1 #Creates points representative of the surface for each subsection

            #Grabs 3 points to make a triangle. 2 from the first foil in the subsection and 1 from the second.
            point1=[first_foil[1][l],first_foil[2][l],first_foil[3][l]]
            point2=[second_foil[1][l],second_foil[2][l],second_foil[3][l]]
            point3=[first_foil[1][l+1],first_foil[2][l+1],first_foil[3][l+1]]
            (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

            #Load information into return dictionary
            push!(meshx,xcoord)
            push!(meshy,ycoord)
            push!(meshz,zcoord)
            push!(mesharea,area)
            push!(meshi,normal[1])
            push!(meshj,normal[2])
            push!(meshk,normal[3])

            #An operation to track the maximum seperation between any two adjacent points
            second_point=[xcoord,ycoord,zcoord]
            if l!=1&&k!=1&&j!=1
              seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
              if seperation>distance
                distance=seperation
              end
            end
            first_point=second_point

            #Sets points on the top of the foil (I think) as having solar panels and those on the bottom without
            if l<85
              push!(meshpaneled,true)
            else
              push!(meshpaneled,false)
            end

            #For symetrical wings create a symmetrical point on the other side of the xz plane
            if aircraft.components[i].symetrical==true
              push!(meshx,xcoord)
              push!(meshy,-ycoord)
              push!(meshz,zcoord)
              push!(mesharea,area)
              push!(meshi,normal[1])
              push!(meshj,-normal[2])
              push!(meshk,normal[3])
              #Sets points on the top of the foil (I think) as having solar panels and those on the bottom without
              if l<85
                push!(meshpaneled,true)
              else
                push!(meshpaneled,false)
              end
            end

            #Cycle points. Keep second ponit from first foil and grab another point from the second foil
            point1=point3
            point3=[second_foil[1][l+1],second_foil[2][l+1],second_foil[3][l+1]]
            (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

            #Load information into return dictionary
            push!(meshx,xcoord)
            push!(meshy,ycoord)
            push!(meshz,zcoord)
            push!(mesharea,area)
            push!(meshi,normal[1])
            push!(meshj,normal[2])
            push!(meshk,normal[3])

            #An operation to track the maximum seperation between any two adjacent points
            second_point=[xcoord,ycoord,zcoord]
            seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
            if seperation>distance
              distance=seperation
            end
            first_point=second_point

            #Sets points on the top of the foil (I think) as having solar panels and those on the bottom without
            if l<85
              push!(meshpaneled,true)
            else
              push!(meshpaneled,false)
            end

            #For symetrical wings create a symmetrical point on the other side of the xz plane
            if aircraft.components[i].symetrical==true
              push!(meshx,xcoord)
              push!(meshy,-ycoord)
              push!(meshz,zcoord)
              push!(mesharea,area)
              push!(meshi,normal[1])
              push!(meshj,-normal[2])
              push!(meshk,normal[3])
              #Sets points on the top of the foil (I think) as having solar panels and those on the bottom without
              if l<85
                push!(meshpaneled,true)
              else
                push!(meshpaneled,false)
              end
            end
          end
          first_foil=second_foil#Sets the second subsection foil as the first subsection foil
        end
        start_foil=end_foil#Sets the second section foil as the first section foil
      end

    #Protocol for components of type fuselage
    elseif aircraft.components[i].is=="fuselage"
      #Intitializes coordinate arrays
      x=Array{Float32}(0)
      y=Array{Float32}(0)
      z=Array{Float32}(0)

      #Create first circular cross section of the fuselage
      for j=0:1:360
        push!(x,aircraft.components[i].sectioncentersx[1])
        push!(y,aircraft.components[i].sectioncentersy[1]+sind(j)*aircraft.components[i].sectionradii[1])
        push!(z,aircraft.components[i].sectioncentersz[1]+cosd(j)*aircraft.components[i].sectionradii[1])
      end
      start_section=[x,y,z]

      n=round(aircraft.components[i].noselength/max_seperation+.5)
      first_section=start_section#Sets first subsection
      for j=1:n-1#Calculate nosecone sections (tip calculated sperately to avoid convergence issues)
        #Intitializes coordinate arrays
        x=Array{Float32}(0)
        y=Array{Float32}(0)
        z=Array{Float32}(0)

        for k=0:1:360
          push!(x,aircraft.components[i].sectioncentersx[1]-sind(j/n*90)*aircraft.components[i].noselength)
          push!(y,aircraft.components[i].sectioncentersy[1]+sind(k)*cosd(j/n*90)*aircraft.components[i].sectionradii[1])
          push!(z,aircraft.components[i].sectioncentersz[1]+cosd(k)*cosd(j/n*90)*aircraft.components[i].sectionradii[1])
        end
        second_section=[x,y,z]

        for k=1:length(second_section[1])-1
          #Grabs 3 points to make a triangle. 2 from the first foil in the subsection and 1 from the second.
          point1=[first_section[1][k],first_section[2][k],first_section[3][k]]
          point2=[second_section[1][k],second_section[2][k],second_section[3][k]]
          point3=[first_section[1][k+1],first_section[2][k+1],first_section[3][k+1]]
          (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

          #Load information into return dictionary
          push!(meshx,xcoord)
          push!(meshy,ycoord)
          push!(meshz,zcoord)
          push!(mesharea,area)
          push!(meshi,normal[1])
          push!(meshj,normal[2])
          push!(meshk,normal[3])
          push!(meshpaneled,false)

          #An operation to track the maximum seperation between any two adjacent points
          second_point=[xcoord,ycoord,zcoord]
          if k!=1&&j!=1
            seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
            if seperation>distance
              distance=seperation
            end
          end
          first_point=second_point

          #Cycle points. Keep second ponit from first foil and grab another point from the second foil
          point1=point3
          point3=[second_section[1][k+1],second_section[2][k+1],second_section[3][k+1]]
          (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

          #Load information into return dictionary
          push!(meshx,xcoord)
          push!(meshy,ycoord)
          push!(meshz,zcoord)
          push!(mesharea,area)
          push!(meshi,normal[1])
          push!(meshj,normal[2])
          push!(meshk,normal[3])
          push!(meshpaneled,false)

          #An operation to track the maximum seperation between any two adjacent points
          second_point=[xcoord,ycoord,zcoord]
          seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
          if seperation>distance
            distance=seperation
          end
          first_point=second_point

        end
        first_section=second_section
      end

      #Perform calculations for end of nosecone
      x=aircraft.components[i].sectioncentersx[1]-aircraft.components[i].noselength
      y=aircraft.components[i].sectioncentersy[1]
      z=aircraft.components[i].sectioncentersz[1]
      tip=[x,y,z]
      for j=1:length(first_section[1])-1
        #Grabs 3 points to make a triangle. 2 from the first foil in the subsection and 1 from the second.
        point1=[first_section[1][j],first_section[2][j],first_section[3][j]]
        #point2=tip
        point3=[first_section[1][j+1],first_section[2][j+1],first_section[3][j+1]]
        (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,tip,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

        #Load information into return dictionary
        push!(meshx,xcoord)
        push!(meshy,ycoord)
        push!(meshz,zcoord)
        push!(mesharea,area)
        push!(meshi,normal[1])
        push!(meshj,normal[2])
        push!(meshk,normal[3])
        push!(meshpaneled,false)
      end

      for j=2:length(aircraft.components[i].sectioncentersx)#Cycles through defined sections
        first_section=start_section#Sets first subsection

        #Intitializes coordinate arrays
        x=Array{Float32}(0)
        y=Array{Float32}(0)
        z=Array{Float32}(0)

        #Create last circular cross section of the section
        for k=0:1:360
          push!(x,aircraft.components[i].sectioncentersx[j])
          push!(y,aircraft.components[i].sectioncentersy[j]+sind(k)*aircraft.components[i].sectionradii[j])
          push!(z,aircraft.components[i].sectioncentersz[j]+cosd(k)*aircraft.components[i].sectionradii[j])
        end
        end_section=[x,y,z]

        n=round((aircraft.components[i].sectioncentersx[j]-aircraft.components[i].sectioncentersx[j-1])/max_seperation+.5)
        for k=1:n
          #Intitializes coordinate arrays
          x=Array{Float32}(0)
          y=Array{Float32}(0)
          z=Array{Float32}(0)

          #Interpolate second subsection cross section using the start and end cross sections of the section... apologies for the redundancy.
          for l=1:length(end_section[1])
            push!(x,start_section[1][l]+k/n*(end_section[1][l]-start_section[1][l]))
            push!(y,start_section[2][l]+k/n*(end_section[2][l]-start_section[2][l]))
            push!(z,start_section[3][l]+k/n*(end_section[3][l]-start_section[3][l]))
          end
          second_section=[x,y,z]

          for l=1:length(end_section[1])-1
            #Grabs 3 points to make a triangle. 2 from the first foil in the subsection and 1 from the second.
            point1=[first_section[1][l],first_section[2][l],first_section[3][l]]
            point2=[second_section[1][l],second_section[2][l],second_section[3][l]]
            point3=[first_section[1][l+1],first_section[2][l+1],first_section[3][l+1]]
            (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

            #Load information into return dictionary
            push!(meshx,xcoord)
            push!(meshy,ycoord)
            push!(meshz,zcoord)
            push!(mesharea,area)
            push!(meshi,normal[1])
            push!(meshj,normal[2])
            push!(meshk,normal[3])
            push!(meshpaneled,false)

            #An operation to track the maximum seperation between any two adjacent points
            second_point=[xcoord,ycoord,zcoord]
            if l!=1&&k!=1&&j!=1
              seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
              if seperation>distance
                distance=seperation
              end
            end
            first_point=second_point

            #Cycle points. Keep second ponit from first foil and grab another point from the second foil
            point1=point3
            point3=[second_section[1][l+1],second_section[2][l+1],second_section[3][l+1]]
            (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

            #Load information into return dictionary
            push!(meshx,xcoord)
            push!(meshy,ycoord)
            push!(meshz,zcoord)
            push!(mesharea,area)
            push!(meshi,normal[1])
            push!(meshj,normal[2])
            push!(meshk,normal[3])
            push!(meshpaneled,false)

            #An operation to track the maximum seperation between any two adjacent points
            second_point=[xcoord,ycoord,zcoord]
            seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
            if seperation>distance
              distance=seperation
            end
            first_point=second_point

          end
          first_section=second_section#Sets the second subsection cross section as the first
        end
        start_section=end_section#Sets the last section cross section as the first
      end

      #Generate hemisphere at end of fuselage
      n=round(aircraft.components[i].sectionradii[length(aircraft.components[i].sectionradii)]/max_seperation+.5)
      for j=1:n-1#Calculate nosecone sections (tip calculated sperately to avoid convergence issues)
        #Intitializes coordinate arrays
        x=Array{Float32}(0)
        y=Array{Float32}(0)
        z=Array{Float32}(0)

        for k=0:1:360
          push!(x,aircraft.components[i].sectioncentersx[length(aircraft.components[i].sectioncentersx)]+sind(j/n*90)*aircraft.components[i].sectionradii[length(aircraft.components[i].sectioncentersx)])
          push!(y,aircraft.components[i].sectioncentersy[length(aircraft.components[i].sectioncentersx)]+sind(k)*cosd(j/n*90)*aircraft.components[i].sectionradii[length(aircraft.components[i].sectioncentersx)])
          push!(z,aircraft.components[i].sectioncentersz[length(aircraft.components[i].sectioncentersx)]+cosd(k)*cosd(j/n*90)*aircraft.components[i].sectionradii[length(aircraft.components[i].sectioncentersx)])
        end
        second_section=[x,y,z]

        for k=1:length(second_section[1])-1
          #Grabs 3 points to make a triangle. 2 from the first foil in the subsection and 1 from the second.
          point1=[first_section[1][k],first_section[2][k],first_section[3][k]]
          point2=[second_section[1][k],second_section[2][k],second_section[3][k]]
          point3=[first_section[1][k+1],first_section[2][k+1],first_section[3][k+1]]
          (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

          #Load information into return dictionary
          push!(meshx,xcoord)
          push!(meshy,ycoord)
          push!(meshz,zcoord)
          push!(mesharea,area)
          push!(meshi,normal[1])
          push!(meshj,normal[2])
          push!(meshk,normal[3])
          push!(meshpaneled,false)

          #An operation to track the maximum seperation between any two adjacent points
          second_point=[xcoord,ycoord,zcoord]
          if k!=1&&j!=1
            seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
            if seperation>distance
              distance=seperation
            end
          end
          first_point=second_point

          #Cycle points. Keep second ponit from first foil and grab another point from the second foil
          point1=point3
          point3=[second_section[1][k+1],second_section[2][k+1],second_section[3][k+1]]
          (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,point2,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

          #Load information into return dictionary
          push!(meshx,xcoord)
          push!(meshy,ycoord)
          push!(meshz,zcoord)
          push!(mesharea,area)
          push!(meshi,normal[1])
          push!(meshj,normal[2])
          push!(meshk,normal[3])
          push!(meshpaneled,false)

          #An operation to track the maximum seperation between any two adjacent points
          second_point=[xcoord,ycoord,zcoord]
          seperation=sqrt((second_point[1]-first_point[1])^2+(second_point[2]-first_point[2])^2+(second_point[3]-first_point[3])^2)
          if seperation>distance
            distance=seperation
          end
          first_point=second_point

        end
        first_section=second_section
      end

      #Perform calculations for end of nosecone
      x=aircraft.components[i].sectioncentersx[length(aircraft.components[i].sectioncentersx)]-aircraft.components[i].sectionradii[length(aircraft.components[i].sectioncentersx)]
      y=aircraft.components[i].sectioncentersy[length(aircraft.components[i].sectioncentersx)]
      z=aircraft.components[i].sectioncentersz[length(aircraft.components[i].sectioncentersx)]
      tip=[x,y,z]
      for j=1:length(first_section[1])-1
        #Grabs 3 points to make a triangle. 2 from the first foil in the subsection and 1 from the second.
        point1=[first_section[1][j],first_section[2][j],first_section[3][j]]
        #point2=tip
        point3=[first_section[1][j+1],first_section[2][j+1],first_section[3][j+1]]
        (xcoord,ycoord,zcoord,area,normal)=triangle_data(point1,tip,point3)#Uses points to create a triangle and calculate the centerpoint, area, and normal vector of the triangle

        #Load information into return dictionary
        push!(meshx,xcoord)
        push!(meshy,ycoord)
        push!(meshz,zcoord)
        push!(mesharea,area)
        push!(meshi,normal[1])
        push!(meshj,normal[2])
        push!(meshk,normal[3])
        push!(meshpaneled,false)
      end
    end
  end
  maxarea=maximum(mesharea)*10^4#Finds the maximum area represented by a single point
  print("Largest Area=$maxarea cm^2\n")#Displays maximum area and maximum distance
  return (surfacemesh(meshx,meshy,meshz,meshi,meshj,meshk,mesharea,meshpaneled,distance,trues(length(meshx))))#Returns dictionary
end

function foil_generator(foil,chord,dihedral,twist,LEx,LEy,LEz)#A function that generates a foil cross section based on the information provided
  foil_data = readdlm("/Users/jcste/OneDrive/Documents/FLOW Lab/Julia Solar Model/"*foil*".dat")#Read in foil geometry
  scaled_foil_data =foil_data*chord#Scale the data

  #Under current definitions a vertical section has a dihedral of 90 degrees. If that changes the following commented out sections will need to be uncommented.
  #=if vertical==true
      cdi = cos(pi-dihedral)
      sdi = sin(pi-diehedral)
  else=#
      cdi = cos(dihedral)
      sdi = sin(dihedral)
  #end
  ctw = cos(twist)
  stw = sin(twist)

  #Rotate according to dihedral
  dihedral_rotation =
    [1  0   0;
     0 cdi sdi;
     0 -sdi cdi ]
  dihedral_foil = ([scaled_foil_data[:,1] zeros(size(scaled_foil_data[:,1])) scaled_foil_data[:,2]]*dihedral_rotation);

  #Rotate according to twist
  twist_rotation =
    [ctw 0 -stw;
     0   1  0;
     stw 0  ctw ]
  twist_foil = (dihedral_foil*twist_rotation)

  #Position foil apropriately
  final_foil=(twist_foil[:,1]+LEx,twist_foil[:,2]+LEy,twist_foil[:,3]+LEz)
  return(final_foil)
end

function triangle_data(point1,point2,point3)#Calculates center coordinates,area, and normal vector of a triangle defined by three points
  #Center point calculated by average x y and z coordinates
  x=(point1[1]+point2[1]+point3[1])/3
  y=(point1[2]+point2[2]+point3[2])/3
  z=(point1[3]+point2[3]+point3[3])/3

  #Area calculated by dot product identity
  ab=sqrt((point2[1]-point1[1])^2+(point2[2]-point1[2])^2+(point2[3]-point1[3])^2)
  ac=sqrt((point3[1]-point1[1])^2+(point3[2]-point1[2])^2+(point3[3]-point1[3])^2)
  AB=[(point2[1]-point1[1]),(point2[2]-point1[2]),(point2[3]-point1[3])]
  AC=[(point3[1]-point1[1]),(point3[2]-point1[2]),(point3[3]-point1[3])]
  θ=cos(dot(AB,AC)/(ab*ac))
  area=.5(ab*ac)sin(θ)

  #Normal vector calculated by cross product
  normal=cross(AB,AC)
  unit_normal=normal/sqrt(normal[1]^2+normal[2]^2+normal[3]^2)

  return(x,y,z,area,unit_normal)
end

using PyCall
using PyPlot
function visualize(airplane) #Simple visualizer function
  Myblue = .85*[.25, .75, 1]
  paneledx=Array{Float32}(0)
  paneledy=Array{Float32}(0)
  paneledz=Array{Float32}(0)
  unpaneledx=Array{Float32}(0)
  unpaneledy=Array{Float32}(0)
  unpaneledz=Array{Float32}(0)
  shadedx=Array{Float32}(0)
  shadedy=Array{Float32}(0)
  shadedz=Array{Float32}(0)

  for i=1:length(airplane.x)
      if airplane.paneled[i]==true&&airplane.lit[i]==true
          push!(paneledx,airplane.x[i])
          push!(paneledy,airplane.y[i])
          push!(paneledz,airplane.z[i])
      elseif airplane[:lit][i]==false
          push!(shadedx,airplane.x[i])
          push!(shadedy,airplane.y[i])
          push!(shadedz,airplane.z[i])
      else
          push!(unpaneledx,airplane.x[i])
          push!(unpaneledy,airplane.y[i])
          push!(unpaneledz,airplane.z[i])
      end
  end
  plot3D(shadedx,shadedy,shadedz,",",markersize=1,color="black")
  plot3D(unpaneledx,unpaneledy,unpaneledz,",",markersize=1,color="gray" )
  plot3D(paneledx,paneledy,paneledz,",",markersize=1,color=Myblue )
  #plot3D(airplane[:x],airplane[:y],airplane[:z],",",markersize=.0001,color=Myblue )
  limit=maximum(airplane.y)*.575#.575 is an experimentally determined scalar meant to increase figure size.
  ax=gca()
  ax[:set_ylim]([-limit,limit])
  ax[:set_zlim]([-limit,limit])
  ax[:set_xlim]([-limit,limit])
  ax[:tick_params]("both",labelsize=0)
  ax[:set_axis_off]()
end
