function generate_gross_cross_section(t, d, b, h, r_inside)

    L = [d, b, h, b, d]
    θ = [π/2, π, 3π/2, 0.0, π/2]
    n = [4, 4, 4, 4, 4]
    n_radius = [3, 3, 3, 3]
    r = r_inside * ones(Float64, length(n_radius))  .+ t  #outside 

    cross_section = CrossSection.generate_open(L, θ, r, n, n_radius)
    unit_node_normals = CrossSection.Tools.calculate_cross_section_unit_node_normals(cross_section)
    x_center, y_center = CrossSection.Tools.get_coords_along_node_normals(cross_section[:,1], cross_section[:, 2], unit_node_normals, -t/2)

    x_center = x_center .- minimum(x_center) .+ t/2
    y_center = y_center .- minimum(y_center) .+ t/2

    coords = (x=x_center, y=y_center)

    return coords

end

function generate_net_cross_section(t, d, b, h, r_inside, h_hole, d_hole, r_hole_stiffener_inside)


    web_depth = h/2 - h_hole/2

    L = [d, b, web_depth, d_hole]
    θ = [π/2, π, 3π/2, 0.0]
    n_radius = [3, 3, 3]
    n = [4, 4, 4, 4]
    r = [r_inside+t, r_inside+t, r_hole_stiffener_inside+t]

    #top half
    cross_section = CrossSection.generate_open(L, θ, r, n, n_radius)
    unit_node_normals = CrossSection.Tools.calculate_cross_section_unit_node_normals(cross_section)
    x_center, y_center = CrossSection.Tools.get_coords_along_node_normals(cross_section[:,1], cross_section[:, 2], unit_node_normals, -t/2)

    x_center = x_center .- minimum(x_center) .+ t/2
    y_center = y_center .- maximum(y_center) .+ h .- t/2

    #reflect

    num_nodes = size(x_center)[1]
    cross_section_bottom = [x_center y_center zeros(Float64, num_nodes)]
    cross_section_bottom = LinesCurvesNodes.rotate_nodes(cross_section_bottom, rotation_axis="x", rotation_center = [t/2, h/2, 0.0], θ=-π)

    coords = (x=[x_center; cross_section_bottom[:,1]], y=[y_center; cross_section_bottom[:,2]])

    return coords

end


function calculate_gross_section_properties(xy_coords)

    coord = [xy_coords.x xy_coords.y]
    num_nodes = size(coord)[1]
    ends = [1:num_nodes-1 2:num_nodes t*ones(Float64, num_nodes-1)] 
    section_props = CUFSM.cutwp_prop2(coord, ends)

    return section_props

end

function calculate_net_section_properties(xy_coords)

    coord = [xy_coords.x xy_coords.y]
    num_nodes = size(coord)[1]
    
    mid_index = Int(num_nodes/2)
    ends_top = [1:mid_index-1 2:mid_index t*ones(Float64, mid_index-1)]
    ends_bottom = [mid_index+1:num_nodes-1 mid_index+2:num_nodes  t*ones(Float64, mid_index-1)]
    ends = [ends_top; ends_bottom]

    section_props = CUFSM.cutwp_prop2(coord, ends)

    return section_props

end


function calculate_Mcrℓ(coords, h, t)

    E = 29500.0
    ν = 0.30

    P = 0.0
    Mxx = 1.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0

    x_center = coords.x
    y_center = coords.y

    lengths = range(1.0, 1.5*h/2, 10)

    curve, shapes, node, elem = CUFSM.open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22)

    half_wavelength = [curve[i][1] for i in eachindex(curve)]
    load_factor = [curve[i][2] for i in eachindex(curve)]
    mode_index = argmin(load_factor)

    Mcrℓ = load_factor[mode_index]

    return Mcrℓ, half_wavelength, load_factor, curve, shapes, node, elem, mode_index

end



function calculate_Mcrd(coords, b, d, h, t)

    CorZ = 0
    bc = b
    dc = d
    θc = 90.0
    μ = 0.30
    ho = h

    Lm = 99999999999.0

    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof=S100AISI.v16.table23131(CorZ,t,bc,dc,θc)
    Lcrd, L = S100AISI.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

    E = 29500.0
    ν = 0.30

    P = 0.0
    Mxx = 1.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0

    x_center = coords.x
    y_center = coords.y

    lengths = range(0.5*Lcrd, 1.5*Lcrd, 10)

    curve, shapes, node, elem = CUFSM.open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22)

    half_wavelength = [curve[i][1] for i in eachindex(curve)]
    load_factor = [curve[i][2] for i in eachindex(curve)]
    mode_index = argmin(load_factor)

    Mcrd = load_factor[mode_index]

    return Mcrd, half_wavelength, load_factor, curve, shapes, node, elem, mode_index

end