using PyCall
using PyPlot
function shade(Model, Conditions,PAR)
  results=Matrix(length(Conditions[:flux]),5)
  for h=1:length(Conditions[:flux])
    close()
    for i=1:length(Model[:name])
      if Model[Model[:name][i]][:type]=="wing"

        geometry=Dict()
        geometry[:LE]=Dict()
        geometry[:TE]=Dict()
        geometry[:CP]=Dict()
        geometry[:LE][:x]=copy(Model[Model[:name][i]][:LE][:x])
        geometry[:LE][:y]=copy(Model[Model[:name][i]][:LE][:y])
        geometry[:LE][:z]=copy(Model[Model[:name][i]][:LE][:z])
        geometry[:TE][:x]=copy(Model[Model[:name][i]][:TE][:x])
        geometry[:TE][:y]=copy(Model[Model[:name][i]][:TE][:y])
        geometry[:TE][:z]=copy(Model[Model[:name][i]][:TE][:z])
        geometry[:CP][:dihedral]=copy(Model[Model[:name][i]][:CP][:dihedral])

        for j=2:length(Model[Model[:name][i]][:LE][:x])
          unshift!(geometry[:LE][:x],Model[Model[:name][i]][:LE][:x][j])
          unshift!(geometry[:LE][:y],-Model[Model[:name][i]][:LE][:y][j])
          unshift!(geometry[:LE][:z],Model[Model[:name][i]][:LE][:z][j])
          unshift!(geometry[:TE][:x],Model[Model[:name][i]][:TE][:x][j])
          unshift!(geometry[:TE][:y],-Model[Model[:name][i]][:TE][:y][j])
          unshift!(geometry[:TE][:z],Model[Model[:name][i]][:TE][:z][j])
          unshift!(geometry[:CP][:dihedral],-Model[Model[:name][i]][:CP][:dihedral][j-1])
        end


        geometry[:CP][:area]=zeros(length(geometry[:LE][:x])-1)
        for j=1:length(geometry[:LE][:x])-1#Calcualte the surface area of each section

          base1=[geometry[:LE][:x][j],geometry[:LE][:y][j],geometry[:LE][:z][j]]
          base2=[geometry[:TE][:x][j],geometry[:TE][:y][j],geometry[:TE][:z][j]]
          top1=[geometry[:LE][:x][j+1],geometry[:LE][:y][j+1],geometry[:LE][:z][j+1]]
          top2=[geometry[:TE][:x][j+1],geometry[:TE][:y][j+1],geometry[:TE][:z][j+1]]

          area=trapezoid_area(base1,base2,top1,top2)
          geometry[:CP][:area][j]=area
        end

        #Load all input angles
        zenith=Conditions[:zenith][h]
        azimuth=Conditions[:azimuth][h]+180
        ψ=-Conditions[:ψ][h]*180/pi
        α=Conditions[:α][h]*180/pi
        ϕ=-Conditions[:ϕ][h]*180/pi

        sun_orientation=-[cosd(zenith)cosd(azimuth) cosd(zenith)sind(azimuth) sind(zenith)]#direction vector pointing to sun

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

        sun_vector=sun_orientation*Rz*Ry*Rx #Rotates sun direction vector about z y and x axis to give direction of sun relative to the plane

        geometry[:CP][:shaded_area]=zeros(length(geometry[:LE][:x])-1)
        geometry[:CP][:shadow_caster]=zeros(length(geometry[:LE][:x])-1)#Boolean stating whether object casts a shadow
        geometry[:CP][:target_section]=zeros(length(geometry[:LE][:x])-1)#Boolean stating whether object is capable of recieving sun light
        geometry[:CP][:percent_perpendicular]=zeros(length(geometry[:LE][:x])-1)
        for j=1:length(geometry[:LE][:x])-1#Loop that cycles through each section
          TE_tip_coord=[geometry[:TE][:x][j+1],geometry[:TE][:y][j+1],geometry[:TE][:z][j+1]]
          LE_base_coord=[geometry[:LE][:x][j],geometry[:LE][:y][j],geometry[:LE][:z][j]]
          TE_base_coord=[geometry[:TE][:x][j],geometry[:TE][:y][j],geometry[:TE][:z][j]]

          TE_vector=unit_vector(TE_tip_coord,TE_base_coord)
          base_vector=unit_vector(LE_base_coord,TE_base_coord)

          section_plane=cross(TE_vector,base_vector)
          percent_perpendicular=-dot(sun_vector,section_plane)
          geometry[:CP][:percent_perpendicular][j]=percent_perpendicular

          if percent_perpendicular<0
            geometry[:CP][:shaded_area][j]=geometry[:CP][:area][j]#If so, the section is entered as being totally shaded
            geometry[:CP][:shadow_caster][j]=1
          else
            geometry[:CP][:target_section][j]=1
          end
        end

        for j=1:length(geometry[:LE][:x])-1
          if geometry[:CP][:shadow_caster][j]==1
            #Coordinates for the corners of the winglet
            LE_tip_coord=[geometry[:LE][:x][j+1],geometry[:LE][:y][j+1],geometry[:LE][:z][j+1]]
            TE_tip_coord=[geometry[:TE][:x][j+1],geometry[:TE][:y][j+1],geometry[:TE][:z][j+1]]
            LE_base_coord=[geometry[:LE][:x][j],geometry[:LE][:y][j],geometry[:LE][:z][j]]
            TE_base_coord=[geometry[:TE][:x][j],geometry[:TE][:y][j],geometry[:TE][:z][j]]

            x=[LE_tip_coord[1],TE_tip_coord[1],TE_base_coord[1],LE_base_coord[1],LE_tip_coord[1]]
            y=[LE_tip_coord[2],TE_tip_coord[2],TE_base_coord[2],LE_base_coord[2],LE_tip_coord[2]]
            z=[LE_tip_coord[3],TE_tip_coord[3],TE_base_coord[3],LE_base_coord[3],LE_tip_coord[3]]


              Myblue = .85*[.25, .75, 1]
              plot_wireframe(x,y,z,color="black")
              limit=maximum(geometry[:LE][:y])*.575#.575 is an experimentally determined scalar meant to increase figure size.
              ax=gca()
              ax[:set_ylim]([-limit,limit])
              ax[:set_zlim]([-limit,limit])
              ax[:set_xlim]([-limit*.1,limit*1.1])
              ax[:tick_params]("both",labelsize=0)
              ax[:set_axis_off]()


            #Find vector representation of winglet's edges
            LE_vector=unit_vector(LE_tip_coord,LE_base_coord)
            TE_vector=unit_vector(TE_tip_coord,TE_base_coord)

            #Find equation of shadow planes (planes that divide lit and shaded portions of surfaces they intersect)
            LE_shadowplane=cross(LE_vector,vec(sun_vector))
            TE_shadowplane=cross(TE_vector,vec(sun_vector))

            for k=1:length(geometry[:LE][:x])-1
              if geometry[:CP][:target_section][k]==1
                #Coordinates of section corners
                left_LE_coord=[geometry[:LE][:x][k],geometry[:LE][:y][k],geometry[:LE][:z][k]]
                left_TE_coord=[geometry[:TE][:x][k],geometry[:TE][:y][k],geometry[:TE][:z][k]]
                right_LE_coord=[geometry[:LE][:x][k+1],geometry[:LE][:y][k+1],geometry[:LE][:z][k+1]]
                right_TE_coord=[geometry[:TE][:x][k+1],geometry[:TE][:y][k+1],geometry[:TE][:z][k+1]]

                x=[left_LE_coord[1],left_TE_coord[1],right_TE_coord[1],right_LE_coord[1],left_LE_coord[1]]
                y=[left_LE_coord[2],left_TE_coord[2],right_TE_coord[2],right_LE_coord[2],left_LE_coord[2]]
                z=[left_LE_coord[3],left_TE_coord[3],right_TE_coord[3],right_LE_coord[3],left_LE_coord[3]]

                plot_wireframe(x,y,z,color=Myblue,linewidth=.2)

                #Find vector representation of section sides
                left_side_vector=unit_vector(left_LE_coord,left_TE_coord)
                right_side_vector=unit_vector(right_LE_coord,right_TE_coord)

                #Find where the shadow intersects the sides of the section (assuming the section is infinite in the x direction)
                front_left_intersect=line_plane_intersect(left_LE_coord,left_side_vector,LE_tip_coord,LE_shadowplane)
                front_right_intersect=line_plane_intersect(right_LE_coord,right_side_vector,LE_tip_coord,LE_shadowplane)
                back_left_intersect=line_plane_intersect(left_LE_coord,left_side_vector,TE_tip_coord,TE_shadowplane)
                back_right_intersect=line_plane_intersect(right_LE_coord,right_side_vector,TE_tip_coord,TE_shadowplane)

                #Find equation for shadowlines on wing face
                LE_shadowline=unit_vector(front_left_intersect,front_right_intersect)
                TE_shadowline=unit_vector(back_left_intersect,back_right_intersect)

                #Calculate where the corners of the shadow would be if they were on this section
                front_left_corner=intersect3d(LE_tip_coord,sun_vector,front_left_intersect,LE_shadowline)
                front_right_corner=intersect3d(LE_base_coord,sun_vector,front_left_intersect,LE_shadowline)
                back_left_corner=intersect3d(TE_tip_coord,sun_vector,back_left_intersect,TE_shadowline)
                back_right_corner=intersect3d(TE_base_coord,sun_vector,back_left_intersect,TE_shadowline)

                Θ=-geometry[:CP][:dihedral][k]
                front_left_corner=x_rotate(front_left_corner,left_LE_coord,Θ)
                front_right_corner=x_rotate(front_right_corner,left_LE_coord,Θ)
                back_left_corner=x_rotate(back_left_corner,left_LE_coord,Θ)
                back_right_corner=x_rotate(back_right_corner,left_LE_coord,Θ)

                right_LE_coord=x_rotate(right_LE_coord,left_LE_coord,Θ)
                right_TE_coord=x_rotate(right_TE_coord,left_LE_coord,Θ)

                #Check if shaded region overlaps target section
                if front_right_corner[2]<=left_LE_coord[2]||front_left_corner[2]>=right_LE_coord[2]
                else
                  LE_shadowline=unit_vector(front_left_corner,front_right_corner)
                  TE_shadowline=unit_vector(back_left_corner,back_right_corner)

                  #Find equations for leading and trailing edge
                  LE_vector=unit_vector(left_LE_coord,right_LE_coord)
                  TE_vector=unit_vector(left_TE_coord,right_TE_coord)

                  finish=right_LE_coord[2]
                  start=left_LE_coord[2]

                  n=(finish-start)/(PAR[:ss])#Set for loop incrementer

                  #Find front and back edges of section at start point
                  s=(start-left_LE_coord[2])/LE_vector[2]
                  wing_front=left_LE_coord+s*LE_vector
                  s=(start-left_TE_coord[2])/TE_vector[2]
                  wing_back=left_TE_coord+s*TE_vector

                  #Find front and back of shadow at start point
                  s=(start-front_left_corner[2])/LE_shadowline[2]
                  left_shadow_front=front_left_corner+s*LE_shadowline
                  s=(start-back_left_corner[2])/TE_shadowline[2]
                  left_shadow_back=back_left_corner+s*TE_shadowline

                  #Set true front and back of shadow
                  if left_shadow_front[1]>wing_back[1]
                    left_shadow_front=wing_back
                  elseif left_shadow_front[1]<wing_front[1]
                    left_shadow_front=wing_front
                  end
                  if left_shadow_back[1]<wing_front[1]
                    left_shadow_back=wing_front
                  elseif left_shadow_back[1]>wing_back[1]
                    left_shadow_back=wing_back
                  end
                  if n==0
                    print("Error")
                  else
                  for l=start+n:n:finish

                    #Find front and back edges of section at end of slice
                    s=(l-left_LE_coord[2])/LE_vector[2]
                    wing_front=left_LE_coord+s*LE_vector
                    s=(l-left_TE_coord[2])/TE_vector[2]
                    wing_back=left_TE_coord+s*TE_vector

                    #Find front and back of shadow at end of slice
                    s=(l-front_left_corner[2])/LE_shadowline[2]
                    right_shadow_front=front_left_corner+s*LE_shadowline
                    s=(l-back_left_corner[2])/TE_shadowline[2]
                    right_shadow_back=back_left_corner+s*TE_shadowline

                    #Set true front and back of shadow
                    if right_shadow_front[1]>wing_back[1]
                      right_shadow_front=wing_back
                    elseif right_shadow_front[1]<wing_front[1]
                      right_shadow_front=wing_front
                    end
                    if right_shadow_back[1]<wing_front[1]
                      right_shadow_back=wing_front
                    elseif right_shadow_back[1]>wing_back[1]
                      right_shadow_back=wing_back
                    end

                    if l>=front_left_corner[2]&&l<=front_right_corner[2]
                    if left_shadow_front!=left_shadow_back&&right_shadow_front!=right_shadow_back
                      front_left=x_rotate(left_shadow_front,left_LE_coord,-Θ)
                      back_left=x_rotate(left_shadow_back,left_LE_coord,-Θ)
                      back_right=x_rotate(right_shadow_back,left_LE_coord,-Θ)
                      front_right=x_rotate(right_shadow_front,left_LE_coord,-Θ)

                      x=[front_left[1],back_left[1],back_right[1],front_right[1],front_left[1]]
                      y=[front_left[2],back_left[2],back_right[2],front_right[2],front_left[2]]
                      z=[front_left[3],back_left[3],back_right[3],front_right[3],front_left[3]]

                      plot_wireframe(x,y,z,color="black",linewidth=.5)

                      #=if k==20
                        print("Caster=$j\n")
                      end=#
                    end

                    area=trapezoid_area(left_shadow_front,left_shadow_back,right_shadow_front,right_shadow_back)
                    geometry[:CP][:shaded_area][k]+=area
                    end

                    left_shadow_front=right_shadow_front
                    left_shadow_back=right_shadow_back
                  end
                  end
                end
              end
            end
          end
        end
      end
      area=0
      shaded_area=0
      flux=0
      for j=1:length(geometry[:CP][:area])
        lit_area=geometry[:CP][:area][j]-geometry[:CP][:shaded_area][j]
          if lit_area<0
            lit_area=0
          end
        area+=geometry[:CP][:area][j]
        shaded_area+=geometry[:CP][:shaded_area][j]
        flux+=geometry[:CP][:percent_perpendicular][j]*lit_area*Conditions[:flux][h]
        if geometry[:CP][:shaded_area][j]/geometry[:CP][:area][j]>1
          print(geometry[:CP][:shaded_area][j]/geometry[:CP][:area][j],"   $j \n")
        end
      end
      percent_shade=shaded_area/area*100
      return percent_shade

      results[h,5]=Conditions[:azimuth][h]
      results[h,4]=Conditions[:ψ][h]
      results[h,3]=flux
      results[h,2]=percent_shade
      results[h,1]=Conditions[:time][h]

      #savefig("image$h")
    end
  end
  if PAR[:write]
    writedlm("data.csv", results,",")
  #=for i=1:4:length(results[:,1])
    close()
    plot(results[:,1],results[:,2],color="blue")
    plot([results[i,1],results[i,1]],[-1,101],color="red")
    ax=gca()
    ax[:set_ylim]([0,100])
    ax[:set_xlim]([9,10])
    title("%Shade vs Time")
    savefig("image$i")

  end=#
end

function intersect3d(point1,line1,point2,line2)#A function to find the interesection of any two 3d lines
  a1=point1[1]
  b1=point1[2]
  c1=point1[3]
  x1=line1[1]
  y1=line1[2]
  z1=line1[3]

  a2=point2[1]
  b2=point2[2]
  c2=point2[3]
  x2=line2[1]
  y2=line2[2]
  z2=line2[3]

  t1=(c1-c2+z1*(b2-b1)/y1)/(z2-z1*(y2/y1))
  intersection=point2+t1*line2
  return intersection
end

function trapezoid_area(base1,base2,top1,top2)#A function to find the area of any 4 sided shape
  δx=base1[1]-base2[1]
  δy=base1[2]-base2[2]
  δz=base1[3]-base2[3]
  base=sqrt(δx^2+δy^2+δz^2)

  δx=top1[1]-top2[1]
  δy=top1[2]-top2[2]
  δz=top1[3]-top2[3]
  top=sqrt(δx^2+δy^2+δz^2)

  δy=base1[2]-top1[2]
  δz=base1[3]-top1[3]
  height1=sqrt(δy^2+δz^2)

  δy=base2[2]-top2[2]
  δz=base2[3]-top2[3]
  height2=sqrt(δy^2+δz^2)

  area=(base+top)*(height1+height2)/4
  return area
end

function unit_vector(point1,point2)#A function that returns the unit vector between two points
  δx=point1[1]-point2[1]
  δy=point1[2]-point2[2]
  δz=point1[3]-point2[3]
  magnitude=sqrt(δx^2+δy^2+δz^2)
  vector=[δx,δy,δz]/magnitude
  return vector
end

function line_plane_intersect(line_point,line_eq,plane_point,plane_eq)
  d=dot(plane_point-line_point,plane_eq)/dot(line_eq,plane_eq)
  intersection=d*line_eq+line_point
  return intersection
end

function x_rotate(point, center_of_rotation,angle) #rotates point about the x-axis superimposed on the center of rotation
  x=point[1]
  y=center_of_rotation[2]+(point[2]-center_of_rotation[2])cos(angle)-(point[3]-center_of_rotation[3])sin(angle)
  z=center_of_rotation[3]+(point[3]-center_of_rotation[3])cos(angle)+(point[2]-center_of_rotation[2])sin(angle)
  new_point=[x,y,z]
  return new_point
end
