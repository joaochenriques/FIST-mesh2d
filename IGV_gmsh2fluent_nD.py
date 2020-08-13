#!/usr/bin/env python
#
# FLUENT set periodic boundary conditions
#

from gmsh2fluent_nD import *
if len( sys.argv ) > 1:
  case_name = sys.argv[1]

else:
  case_name = "test1"
#===============================================================================
if __name__ == "__main__":

    # Volume id used in fluent
    export_vol_id = volume( "profiles", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior-domain",                 2, bc_type_interior, 0 ),
        surface( "periodic_bc_up",                  3, bc_type_interface, 32  ),
        surface( "outlet_prof",                     4, bc_type_wall, 2  ),
        surface( "periodic_bc_down",                5, bc_type_interface, 4  ),
        surface( "prof_in",                         6, bc_type_wall, 8  ),
        surface( "blade_bl1",                       7, bc_type_wall, 1  ),
        surface( "blade_bl2",                       8, bc_type_wall, 128  ),
        surface( "blade_bl3",                       9, bc_type_wall, 256  ),
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_int.msh" % case_name, scale = 1000.0 )
    del mesh2d

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# BLADE 1
    # Volume id used in fluent
    export_vol_id = volume( "bl1_blk1", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior_bl11",                      2, bc_type_interior, 0 ),
        surface( "outer_bl11",                         3, bc_type_wall, 1 ),
        surface( "blade_wall_intr",                    4, bc_type_wall, 2  ),     # linha 2 corresponde as superficies, logo tem intra e extradorso!!!!!!!!!!!!!
        surface( "blade_wall_extr",                    5, bc_type_wall, 64  ),
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_BL1_blk1_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_bl11.msh" % case_name, scale = 1000.0 )
    del mesh2d
    
     # Volume id used in fluent
    export_vol_id = volume( "bl1_blk2", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior_bl12",                      2, bc_type_interior, 0 ),
        surface( "outer_bl12",                         3, bc_type_wall, 1 ),    
        surface( "inner_bl12",                         4, bc_type_wall, 2  ),   
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_BL1_blk2_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_bl12.msh" % case_name, scale = 1000.0 )
    del mesh2d    

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# BLADE 2    
    # Volume id used in fluent
    export_vol_id = volume( "bl2_blk1", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior_bl21",                      2, bc_type_interior, 0 ),
        surface( "outer_bl21",                         3, bc_type_wall, 1 ),
        surface( "blade_wall_intr_flap",                  4, bc_type_wall, 2  ),
        surface( "blade_wall_extr_flap",                  5, bc_type_wall, 64  ),
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_BL2_blk1_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_bl21.msh" % case_name, scale = 1000.0 )
    del mesh2d
    
         # Volume id used in fluent
    export_vol_id = volume( "bl2_blk2", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior_bl22",                      2, bc_type_interior, 0 ),
        surface( "outer_bl22",                         3, bc_type_wall, 1 ),
        surface( "inner_bl22",                         4, bc_type_wall, 2  ),
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_BL2_blk2_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_bl22.msh" % case_name, scale = 1000.0 )
    del mesh2d    
    
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # Volume id used in fluent
    export_vol_id = volume( "esteira", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior_esteira",                 2, bc_type_interior, 0 ),
        surface( "in_est",                           3, bc_type_wall, 2  ),
        surface( "bottom_est",                       4, bc_type_interface, 64 ),
        surface( "pressure_outlet",                  5, bc_type_pressure_outlet, 8 ),
        surface( "top_est",                          6, bc_type_interface, 16 ),
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_est_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_est.msh" % case_name, scale = 1000.0 )
    del mesh2d

      # Volume id used in fluent
    export_vol_id = volume( "entrada", 1 )

    # Boundary conditions used in fluent
    export_srfs = (
        surface( "interior_entrada",                 2, bc_type_interior, 0 ),
        surface( "velocity_in",                      3, bc_type_velocity_inlet, 2  ),
        surface( "bottom_entrada",                   4, bc_type_interface, 64 ),
        surface( "entrada_outlet",                   5, bc_type_wall, 8 ),
        surface( "top_entrada",                      6, bc_type_interface, 16 ),
    )

    mesh2d = mesh(2)
    mesh2d.read_gmsh( "./%s_entr_gmsh.msh" % case_name )
    mesh2d.write_fluent( export_vol_id, export_srfs,
                         "./%s_entr.msh" % case_name, scale = 1000.0 )
    del mesh2d
print "END"
    
#=============================================================================

#sys.exit(0)
##==EOF==========================================================================