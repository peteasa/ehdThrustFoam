#####################################################################
# inspired by Ole Petersen from his blog post
# https://peteole.github.io/blog/openfoam_gmsh_2d_aerofoil/index.html
#####################################################################

y(x)= 0.594689181*(0.298222773*sqrt(x) - 0.127125232*x - 0.357907906*x^2 + 0.291984971*x^3 - 0.105174606*x^4)

# calculation based on 0.1m airfoil length
# if changed change the mesh_size to match
fractairfoillength = 0.1
boundarylength = 0.1
n_points_per_side=200
top_elect_start = 0
top_elect_end = 51
top_tail_start = n_points_per_side - 30
top_tail_end = n_points_per_side - top_elect_start + 1
bottom_elect_start = n_points_per_side - top_elect_end + 1
bottom_elect_end = n_points_per_side - top_elect_start + 1
bottom_tail_start = n_points_per_side - top_tail_end + 1
bottom_tail_end = n_points_per_side - top_tail_start + 1
n_thrust_surface_points = 2*n_points_per_side
thrust_surface_length = 0.3
thrust_surface_height = 0.1
bounding_radius=0.10
electrode_position=-0.0037
electrode_radius=0.0007
far_field_mesh_size=0.004
pelement_mesh_size=0.0001
electrode_mesh_size=0.0001
thrust_surface_mesh_size=0.001
run(`./set_config_property.sh name unstructured_gmsh_$(electrode_mesh_size)_$(far_field_mesh_size)`)
run(`./set_config_property.sh separation $(- electrode_position - electrode_radius)`)
touch("ehdthrust.geo")
open("ehdthrust.geo","w")do io
    println(io,"""SetFactory("OpenCASCADE");""")
    println(io,"""meshThickness = 1.0;""")

    for i in 1:n_points_per_side
        x=(i-1)/(n_points_per_side)
        #mesh_size=n_thrust_surface_points#0.05#1*fractairfoillength*x*(1-x)+0.005
        if top_elect_start < i && i < top_elect_end
           println(io,"Point($i) = {$(fractairfoillength*x), $(fractairfoillength*y(x)), 0.0,$electrode_mesh_size};")
        elseif top_tail_start < i && i < top_tail_end
           println(io,"Point($i) = {$(fractairfoillength*x), $(fractairfoillength*y(x)), 0.0,$electrode_mesh_size};")
        else
           println(io,"Point($i) = {$(fractairfoillength*x), $(fractairfoillength*y(x)), 0.0,$thrust_surface_mesh_size};")
        end
    end
    for i in 1:n_points_per_side
        x=1-(i-1)/(n_points_per_side)
        #mesh_size=n_thrust_surface_points#0.05#1*fractairfoillength*x*(1-x)+0.005
        if bottom_elect_start < i && i < bottom_elect_end
           println(io,"Point($(i+n_points_per_side)) = {$(fractairfoillength*x), $(-fractairfoillength*y(x)), 0.0,$electrode_mesh_size};")
        elseif bottom_tail_start < i && i < bottom_tail_end
           println(io,"Point($(i+n_points_per_side)) = {$(fractairfoillength*x), $(-fractairfoillength*y(x)), 0.0,$electrode_mesh_size};")
        else
           println(io,"Point($(i+n_points_per_side)) = {$(fractairfoillength*x), $(-fractairfoillength*y(x)), 0.0,$thrust_surface_mesh_size};")
        end
    end

    println(io,"Point($(n_thrust_surface_points+1)) = {0.0, $(bounding_radius), 0.0, $far_field_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+2)) = {0.0, $(-bounding_radius), 0.0, $far_field_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+3)) = {$(-bounding_radius), 0.0, 0.0, $far_field_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+4)) = {$(bounding_radius+boundarylength*1), $(bounding_radius), 0.0, $far_field_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+5)) = {$(bounding_radius+boundarylength*1), $(-bounding_radius), 0.0, $far_field_mesh_size};")

    # now points for the PELEMENT feature
    println(io,"Point($(n_thrust_surface_points+6)) =  {$(electrode_position - electrode_radius), 0.0, 0.0, $pelement_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+7)) =  {$electrode_position, +$(electrode_radius), 0.0, $pelement_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+8)) =  {$(electrode_position + electrode_radius), 0.0, 0.0, $pelement_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+9)) =  {$electrode_position, -$(electrode_radius), 0.0, $pelement_mesh_size};")
    println(io,"Point($(n_thrust_surface_points+10)) =  {$electrode_position, 0.0, 0.0, $pelement_mesh_size};")

    # number of surfaces created at start of extrude
    n_extrude = 2
    # define the lines that make up the outer boundary Note work in consistent anti clockwise direction
    # the normal to the extruded surfaces will then point inwards
    # in practice paraFoam shows the normals point in the opposite direction
    inlet_upper = n_extrude
    println(io,"Circle($(n_thrust_surface_points+1)) = {$(n_thrust_surface_points+2), 1, $(n_thrust_surface_points+3)};")
    inlet_lower = inlet_upper + 1
    println(io,"Circle($(n_thrust_surface_points+2)) = {$(n_thrust_surface_points+3), 1, $(n_thrust_surface_points+1)};")
    wall_lower = inlet_lower + 1
    println(io,"Line($(n_thrust_surface_points+3)) = {$(n_thrust_surface_points+1), $(n_thrust_surface_points+4)};")
    outlet = wall_lower + 1
    println(io,"Line($(n_thrust_surface_points+4)) = {$(n_thrust_surface_points+4), $(n_thrust_surface_points+5)};")
    wall_upper = outlet + 1
    println(io,"Line($(n_thrust_surface_points+5)) = {$(n_thrust_surface_points+5), $(n_thrust_surface_points+2)};")

    # PELEMENT feature
    println(io,"Circle($(n_thrust_surface_points+6)) = {$(n_thrust_surface_points+6), $(n_thrust_surface_points+10), $(n_thrust_surface_points+7)};")
    println(io,"Circle($(n_thrust_surface_points+7)) = {$(n_thrust_surface_points+7), $(n_thrust_surface_points+10), $(n_thrust_surface_points+8)};")
    println(io,"Circle($(n_thrust_surface_points+8)) = {$(n_thrust_surface_points+8), $(n_thrust_surface_points+10), $(n_thrust_surface_points+9)};")
    println(io,"Circle($(n_thrust_surface_points+9)) = {$(n_thrust_surface_points+9), $(n_thrust_surface_points+10), $(n_thrust_surface_points+6)};")

    c = 2
    # create the outer boundary work round the boundary
    outer = c
    println(io,"Curve Loop($c) = {$(join((n_thrust_surface_points+1):(n_thrust_surface_points+5),","))};") ; c += 1

    l = 1
    # create the airfoil work round the airfoil Note work in consistent clockwise direction
    # the normal to the extruded surfaces will then point away from the airfoil
    # in practice paraFoam shows the normals point in the opposite direction
    upper_nu = wall_upper + 1
    println(io,"Spline($l) = {$(top_elect_start+1):$top_elect_end};") ; l += 1
    upper_a = upper_nu + 1
    println(io,"Spline($l) = {$top_elect_end:$top_tail_start};") ; l += 1
    upper_t = upper_a + 1
    println(io,"Spline($l) = {$top_tail_start:$top_tail_end};") ; l += 1
    lower_t = upper_t + 1
    println(io,"Spline($l) = {$(n_points_per_side+bottom_tail_start+1):$(n_points_per_side+bottom_tail_end)};") ; l += 1
    lower_b = lower_t + 1
    println(io,"Spline($l) = {$(n_points_per_side+bottom_tail_end):$(n_points_per_side+bottom_elect_start)};") ; l += 1
    lower_nu = lower_b + 1
    println(io,"Spline($l) = {$(n_points_per_side+bottom_elect_start):$(2*n_points_per_side), 1};") ; l += 1
    # println(io,"Line($l) = {$n_thrust_surface_points, 1};") ; l += 1

    airfoil = c
    println(io,"Curve Loop($c) = {$(join(1:(l-1),","))};") ; c += 1

    element = c
    pelem = lower_nu + 1 # first of four circle segments
    println(io,"Curve Loop($c) = {$(n_thrust_surface_points+6), $(n_thrust_surface_points+7), $(n_thrust_surface_points+8), $(n_thrust_surface_points+9)};") ; c += 1

    println(io,"Plane Surface(1) = {$outer, $airfoil, $element};")

    println(io,"""
    surfaceVector[] = Extrude {0, 0, meshThickness} {
        Surface{1};
        Layers{1};
        Recombine;
    };
    For s In {0:#surfaceVector[]-1}
        Printf("extruded surface %g: %g ", s, surfaceVector[s]);
    EndFor
    Physical Volume("internalField") = surfaceVector[1];
    Physical Surface("frontAndBackPlanes") = {1, surfaceVector[0]};
    Physical Surface("INLET")={surfaceVector[$inlet_upper],surfaceVector[$inlet_lower]};
    Physical Surface("OUTLET")={surfaceVector[$outlet]};
    Physical Surface("WALL")={surfaceVector[$wall_lower],surfaceVector[$wall_upper]};
    """)
    println(io,"""Physical Surface("NELEMENT")={surfaceVector[$upper_nu], surfaceVector[$upper_a], surfaceVector[$upper_t], surfaceVector[$lower_t], surfaceVector[$lower_b], surfaceVector[$lower_nu]};""")
    println(io,"""Physical Surface("PELEMENT")={surfaceVector[$pelem], surfaceVector[$(pelem+1)], surfaceVector[$(pelem+2)], surfaceVector[$(pelem+3)]};""")
    println(io,"""Recombine Surface{1};""")

    #println(io,"Physical Curve(\"BOUNDARY\") = {$(join((n_thrust_surface_points+1):(n_thrust_surface_points+5),","))};")
end
