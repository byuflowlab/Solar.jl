# ----- Trajectory Types ------ #
"""
   aircrafttrajectory

Holds aircraft orientation information for each time step
"""
struct aircrafttrajectory
    time::Array{Float64,1}
    roll::Array{Float64,1}
    pitch::Array{Float64,1}
    yaw::Array{Float64,1}
end

# ----- Solar Panel Types ----- #
"
    panelgeometry

Describe the geometry/orientation of a solar panel array.
"
struct panelgeometry
    nodes::Array{Array{Float64,1},1}
    centroid::Array{Float64,1}
    area::Float64
    normal::Array{Float64,1}
end

function panelgeometry(area, normal)
    nodes = Array{Array{Float64,1},1}(0)
    centroid = Array{Float64,1}(0)
    panelgeometry(nodes, centroid, area, normal)
end

function panelgeometry(nodes, area, normal)
    x = [nodes[i][1] for i in eachindex(nodes)]
    y = [nodes[i][2] for i in eachindex(nodes)]
    z = [nodes[i][3] for i in eachindex(nodes)]
    centroid = [mean(x),mean(y),mean(z)]
    panelgeometry(nodes, centroid, area, normal)
end



"""
    panelparameters

"""
struct panelparameters
    #TODO talk to Nathaniel about what needs to go in here as far as the power stuff is concerned
end #panelphysicsin

"""
    panelenergy

Define flux, power, and total energy arrays for a solar panel network.
"""
struct panelenergy
    shadow::Array{Bool,2}
    flux::Array{Float64,2}
    power::Array{Float64,2}
    energy::Array{Float64,2}
end #panelphysicsout

function panelenergy(flux, power, energy)
    shadow = zeros(Bool, size(flux)...)
    panelenergy(shadow, flux, power, energy)
end

# ----- Sunshine Types ----- #
"""
    sunshine

Define time, sun angles (azimuth, zenith), and solar flux arrays.
"""
struct sunshine
    time::Array{Float64,1}
    azimuth::Array{Float64,1}
    zenith::Array{Float64,1}
    flux::Array{Float64,1}
end

# ----- Aircraft Geometry Types ----- #
struct winggeometry
    is::String
    LEx::Array{Float64,1}
    LEy::Array{Float64,1}
    LEz::Array{Float64,1}
    dihedral::Array{Float64,1}
    twist::Array{Float64,1}
    chord::Array{Float64,1}
    symetrical::Bool
    vertical::Bool
    foil::Array{String,1}
end
struct fuselagegeometry
    is::String
    noseposex::Float64
    noseposey::Float64
    noseposez::Float64
    noselength::Float64
    sectioncentersx::Array{Float64,1}
    sectioncentersy::Array{Float64,1}
    sectioncentersz::Array{Float64,1}
    sectionradii::Array{Float64,1}
end
struct aircraftgeometry
    components::Array{Any,1}
end
function definegeometry()
    #define fusealge geometry
    is="fuselage"
    npx=-5
    npy=0
    npz=-2
    nl=4
    scx=[-1.5,-.563,1.69,2.69,3.69,4.69,5.69,6.69,7.69,8.69,9.69,10.69]
    scy=[0,0,0,0,0,0,0,0,0,0,0,0]
    scz=[-1,-1,-1,-1,-1,-1,-.9,-.8,-.7,-.6,-.5,-.5]
    r=[1,1,1,1,1,1,.9,.8,.7,.6,.5,.5]
    fuselage=fuselagegeometry(is,npx,npy,npz,nl,scx,scy,scz,r)

    #define horizontal stabilizer geometry
    is="wing"
    lex=[8.69,9.19,9.69]
    ley=[0,2,4]
    lez=[0,0,0]
    d=[0,0]
    t=[0,0,0]
    c=[2,1.5,1]
    s=true
    v=false
    f=["naca2412","naca2412","naca2412"]
    horizontal_stabilizer=winggeometry(is,lex,ley,lez,d,t,c,s,v,f)

    #define vertical stabilizer geometry
    is="wing"
    lex=[8.69,9.69]
    ley=[0,0]
    lez=[0,2]
    d=[pi/2]
    t=[0,0]
    c=[2,1]
    s=false
    v=true
    f=["naca2412","naca2412"]
    vertical_stabilizer=winggeometry(is,lex,ley,lez,d,t,c,s,v,f)

    #define main wing geometry
    is="wing"
    lex=[-0.5625,0.057904013549477606,0.7858444820818824,1.513784950614287,2.2417254191466918,2.9696658876790964,3.6976063562115007,4.425546824743906,5.15348729327631,5.881427761808715,6.609368230341119,7.3373086988735245,8.112047487083235,8.886786275292943]
    ley=[0.00,1.1207190353532137,3.0525706879313503,4.984422340509487,6.916273993087623,8.84812564566576,10.779977298243896,12.711828950822033,14.64368060340017,16.575532255978306,18.50738390855644,20.43923556113458,20.43923556113458,20.43923556113458]
    lez=[0.00,0.09805021059111543,0.6156883007961569,1.1333263910011984,1.6509644812062398,2.1686025714112813,2.686240661616323,3.2038787518213643,3.7215168420264058,4.239154932231448,4.756793022436489,5.274431112641531,6.274431112641531,7.274431112641531]
    d=[0.08726646259971647,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,0.2617993877991494,1.5707963267948966,1.5707963267948966]
    t=[0.00,0.00,-0.07853981633974483,-0.15707963267948966,-0.23561944901923448,-0.3141592653589793,-0.39269908169872414,-0.47123889803846897,-0.5497787143782138,-0.6283185307179586,-0.7068583470577035,-0.7853981633974483,-0.39269908169872414,0.00]
    c=[2.25,1.40625,1.40625,1.40625,1.40625,1.40625,1.40625,1.40625,1.40625,1.40625,1.40625,1.40625,1.108125,0.8099999999999999]
    s=true
    v=false
    f=["sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","sg6043","naca2412"]
    main_wing=winggeometry(is,lex,ley,lez,d,t,c,s,v,f)

    return aircraftgeometry([main_wing,vertical_stabilizer,horizontal_stabilizer,fuselage])
end

struct surfacemesh
    x::Array{Float32,1}
    y::Array{Float32,1}
    z::Array{Float32,1}
    i::Array{Float32,1}
    j::Array{Float32,1}
    k::Array{Float32,1}
    area::Array{Float32,1}
    paneled::Array{Bool}
    precision::Float32
end
