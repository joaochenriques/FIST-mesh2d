#!/usr/bin/env python
#===============================================================================
# This is not a generic code. Suported cell types:
# * 2D - lines / triangles / quadrangles
# * 3D - triangles / tetrahedras
#===============================================================================
import string, re, sys, math

#===============================================================================
# FLUENT BC
bc_type_interior          = (  2, "interior" )
bc_type_wall              = (  3, "wall" )
bc_type_pressure_inlet    = (  4, "pressure-inlet" )
bc_type_pressure_outlet   = (  5, "pressure-outlet" )
bc_type_symmetry          = (  7, "symmetry" )
bc_type_periodic_shadow   = (  8, "shadow" )
bc_type_pressure_farfield = (  9, "pressure-far-field" )
bc_type_velocity_inlet    = ( 10, "velocity-inlet" )
bc_type_periodic          = ( 12, "periodic" )
bc_type_fan               = ( 14, "intake-fan" )
bc_type_mass_flow_inlet   = ( 20, "mass-flow-inlet" )
bc_type_interface         = ( 24, "interface" )
bc_type_outflow           = ( 36, "outflow" )
bc_type_axis              = ( 37, "axis" )

#===============================================================================
# FLUENT BC not supported
# bc_type_hanging_node      = ( 31, "hanging-node" )

#===============================================================================
re_ws         = r"(\s*)"
re_begin_node = r"(\s*)(\$NOD)(\s*)"
re_end_node   = r"(\s*)(\$ENDNOD)(\s*)"
re_begin_elm  = r"(\s*)(\$ELM)(\s*)"
re_end_elm    = r"(\s*)(\$ENDELM)(\s*)"
re_begin_per  = r"(\s*)(\$PERNODES)(\s*)"
re_end_per    = r"(\s*)(\$ENDPERNODES)(\s*)"

#===============================================================================
def print_tokens( tokens ):
  if tokens:
    for i in range( 0, len( tokens.groups() ) ):
      print i, ":", tokens.group(i)

#===============================================================================
def exit_error( code ):
  print "ERROR@" + code
  import traceback
  traceback.print_stack()
  exit(1)

#===============================================================================
def hash_func( lst ):
  n = len( lst )
  lst = list(lst)
  lst.sort()

  if n == 3:   # most frequent cases are coded
    fmt = "%08x-%08x-%08x"
    hash = ( fmt % tuple(lst) )
    return hash

  elif n == 4:
    fmt = "%08x-%08x-%08x-%08x"
    hash = ( fmt % tuple(lst) )
    return hash

  elif n == 2:
    fmt = "%08x-%08x"
    hash = ( fmt % tuple(lst) )
    return hash

  else: # generic case
    fmt = "%08x-"*(n-1) + "%08x"
    hash = ( fmt % tuple(lst) )
    return hash

#===============================================================================
class node:

  def __init__( self, node_number, x, y, z ):
    self.node_number = node_number
    self.point = ( x, y, z )
    self.fl_number = -1          # fluent node number

#===============================================================================
class cell:

  ## since fluent reads the elements by "faces", the face definition depends
  ## on the dimension of the problem
  ##
  @classmethod
  def set_dimension( cls, dimension ):
    assert( dimension == 2 or dimension == 3 )
    cls.__dim = dimension

  ##----------------------------------------------------------------------------
  ## left handed 2D mesh
  TRI_FAC2D   = ( ( 0, 1 ), ( 1, 2 ), ( 2, 0 ) )
  EDGE_FAC2D  = ( 0, 1 )
  QUAD_FAC2D  = ( ( 0, 1 ), ( 1, 2 ), ( 2, 3 ), ( 3, 0 ) )

  ## right handed 2D mesh
  #~ EDGE_FAC2D  = ( 1, 0 )
  #~ TRI_FAC2D   = ( ( 0, 2 ), ( 2, 1 ), ( 1, 0 ) )
  #~ QUAD_FAC2D  = ( ( 0, 3 ), ( 3, 2 ), ( 2, 1 ), ( 1, 0 ) )

  ##----------------------------------------------------------------------------
  TRI_FAC3D = ( ( 0,1,2 ), )

  QUAD_FAC3D = ( ( 0,1,2,3 ), )

  #~ TETRA_FAC3D = ( ( 1,2,3 ), ( 0,3,2 ), ( 0,1,3 ), ( 0,2,1 ) ) # RH
  TETRA_FAC3D = ( ( 1,3,2 ), ( 0,2,3 ), ( 0,3,1 ), ( 0,1,2 ) ) # LH

  #~ PYRAM_FAC3D = ( ( 3,2,1,0 ), ( 0,1,4 ), ( 1,2,4 ), ( 2,3,4 ), ( 3,0,4 ) ) # RH
  PYRAM_FAC3D = ( ( 0,1,2,3 ), ( 4,1,0 ), ( 4,2,1 ), ( 4,3,2 ), ( 4,0,3 ) ) # LH

  #~ PRISM_FAC3D = ( ( 2,1,0 ), ( 3,4,5 ), ( 5,4,1,2 ), ( 0,3,5,2 ), # RH
                  #~ ( 0,1,4,3 ) )
  PRISM_FAC3D = ( ( 0,1,2 ), ( 5,4,3 ), ( 2,1,4,5 ), ( 2,5,3,0 ), # LH
                  ( 3,4,1,0 ) )

  HEX_FAC3D = ( ( 3,2,1,0 ), ( 4,5,6,7 ), ( 2,3,7,6 ), ( 0,1,5,4 ), # RH
                ( 0,4,7,3 ), ( 1,2,6,5 ) )
  #~ HEX_FAC3D = ( ( 0,1,2,3 ), ( 7,6,5,4 ), ( 6,7,3,2 ), ( 4,5,1,0 ), # LH
                #~ ( 3,7,4,0 ), ( 5,6,2,1 ) )

  ##----------------------------------------------------------------------------
  def __hash_func( self, nodes_indx ):
    lst = []
    for n in nodes_indx:
      lst.append( self.nodes_lst[ n ].node_number )
    return hash_func( lst )

  ##----------------------------------------------------------------------------
  def __init__( self, elm_number, elm_type, reg_phys, reg_elem, nodes_lst ):
    self.elm_number = elm_number
    self.elm_type   = elm_type
    self.reg_phys   = reg_phys
    self.reg_elem   = reg_elem
    self.nodes_lst  = nodes_lst
    self.fl_number  = -1          # fluent cell number

    if cell.__dim == 2: # the faces hash is dimension dependent

      if elm_type == 1:
        self.faces_hash = ( self.__hash_func( self.EDGE_FAC2D ),  )

      elif elm_type == 2:
        self.faces_hash = ( self.__hash_func( self.TRI_FAC2D[0] ),
                            self.__hash_func( self.TRI_FAC2D[1] ),
                            self.__hash_func( self.TRI_FAC2D[2] ) )

      elif elm_type == 3:
        self.faces_hash = ( self.__hash_func( self.QUAD_FAC2D[0] ),
                            self.__hash_func( self.QUAD_FAC2D[1] ),
                            self.__hash_func( self.QUAD_FAC2D[2] ),
                            self.__hash_func( self.QUAD_FAC2D[3] ) )

      else:
        exit_error( "invalid cell type 2D" )

    else:

      if elm_type == 2:
        self.faces_hash = ( self.__hash_func( self.TRI_FAC3D[0] ), )

      elif elm_type == 3:
        self.faces_hash = ( self.__hash_func( self.QUAD_FAC3D[0] ), )

      elif elm_type == 4:
        self.faces_hash = ( self.__hash_func( self.TETRA_FAC3D[0] ),
                            self.__hash_func( self.TETRA_FAC3D[1] ),
                            self.__hash_func( self.TETRA_FAC3D[2] ),
                            self.__hash_func( self.TETRA_FAC3D[3] ) )

      elif elm_type == 5:
        self.faces_hash = ( self.__hash_func( self.HEX_FAC3D[0] ),
                            self.__hash_func( self.HEX_FAC3D[1] ),
                            self.__hash_func( self.HEX_FAC3D[2] ),
                            self.__hash_func( self.HEX_FAC3D[3] ),
                            self.__hash_func( self.HEX_FAC3D[4] ),
                            self.__hash_func( self.HEX_FAC3D[5] ) )

      elif elm_type == 6:
        self.faces_hash = ( self.__hash_func( self.PRISM_FAC3D[0] ),
                            self.__hash_func( self.PRISM_FAC3D[1] ),
                            self.__hash_func( self.PRISM_FAC3D[2] ),
                            self.__hash_func( self.PRISM_FAC3D[3] ),
                            self.__hash_func( self.PRISM_FAC3D[4] ) )

      elif elm_type == 7:
        self.faces_hash = ( self.__hash_func( self.PYRAM_FAC3D[0] ),
                            self.__hash_func( self.PYRAM_FAC3D[1] ),
                            self.__hash_func( self.PYRAM_FAC3D[2] ),
                            self.__hash_func( self.PYRAM_FAC3D[3] ),
                            self.__hash_func( self.PYRAM_FAC3D[4] ) )

      else:
        print elm_type,
        exit_error( "invalid cell type 3D" )

  ##----------------------------------------------------------------------------
  def face_hash( self, f ):
    return self.faces_hash[f]

  ##----------------------------------------------------------------------------
  def face_nodes( self, f ):

    if cell.__dim == 2: # the faces nodes is dimension dependent

      if self.elm_type == 1:
        return tuple( self.nodes_lst )

      elif self.elm_type == 2:
        return ( self.nodes_lst[ self.TRI_FAC2D[f][0] ],
                 self.nodes_lst[ self.TRI_FAC2D[f][1] ] )

      elif self.elm_type == 3:
        return ( self.nodes_lst[ self.QUAD_FAC2D[f][0] ],
                 self.nodes_lst[ self.QUAD_FAC2D[f][1] ] )

      else:
        exit_error( "invalid cell type 2D" )

    else:

      if self.elm_type == 2:
        return tuple( self.nodes_lst )

      elif self.elm_type == 3:
        return ( self.nodes_lst )

      elif self.elm_type == 4:
        return ( self.nodes_lst[ self.TETRA_FAC3D[f][0] ],
                 self.nodes_lst[ self.TETRA_FAC3D[f][1] ],
                 self.nodes_lst[ self.TETRA_FAC3D[f][2] ] )

      elif self.elm_type == 5:
        return ( self.nodes_lst[ self.HEX_FAC3D[f][0] ],
                 self.nodes_lst[ self.HEX_FAC3D[f][1] ],
                 self.nodes_lst[ self.HEX_FAC3D[f][2] ],
                 self.nodes_lst[ self.HEX_FAC3D[f][3] ] )

      elif self.elm_type == 6:
        indx_lst = self.PRISM_FAC3D[f]
        nodes_lst = []
        for i in indx_lst:
          nodes_lst.append( self.nodes_lst[i] )
        return tuple( nodes_lst )

      elif self.elm_type == 7:
        indx_lst = self.PYRAM_FAC3D[f]
        nodes_lst = []
        for i in indx_lst:
          nodes_lst.append( self.nodes_lst[i] )
        return tuple( nodes_lst )

      else:
        exit_error( "invalid cell type 3D" )

  ##----------------------------------------------------------------------------
  def num_faces( self ):
    return len( self.faces_hash )

#===============================================================================
class surface:

  def __init__( self, fl_name, fl_id, fl_bc_type, gmsh_bc_type ):
    self.fl_name = fl_name
    self.fl_id = fl_id
    self.fl_bc_type = fl_bc_type
    self.gmsh_bc_type = gmsh_bc_type

#===============================================================================
class volume:

  def __init__( self, fl_name, fl_id ):
    self.fl_name = fl_name
    self.fl_id = fl_id

#===============================================================================
class face_cnx:
  def __init__( self, right_cell, right_face ):

    self.right_face = right_face
    self.right_cell = right_cell  # face number

    self.left_cell = None
    self.left_face = None         # face number

    self.bc_type = None

#===============================================================================
class mesh:

  ##----------------------------------------------------------------------------
  def __init__( self, dim ):
    self.mesh_dim = dim
    cell.set_dimension( dim )

    self.nodes_map = {}
    self.nodes_lst = [] # to be sorted

    self.cells_map = {}
    self.cells_lst = [] # to be sorted

    self.file_lines = None
    self.file_lineno = None

    self.periodic_map = {}
    self.cell_type = {}

  ##----------------------------------------------------------------------------
  def __gmsh_read_nodes( self ):
    rgx = re.compile( re_begin_node )
    tokens = rgx.search( self.file_lines[ self.file_lineno ] )
    if tokens:
      self.file_lineno += 1
      fields = string.split( self.file_lines[ self.file_lineno ] )
      num_nodes = int( fields[0] )
      self.file_lineno += 1
      for j in range( 0, num_nodes ):
        fields = string.split( self.file_lines[ self.file_lineno ] )
        node_number = int( fields[0] )
        x = float( fields[1] )
        y = float( fields[2] )
        z = float( fields[3] )
        new_node = node( node_number, x, y, z )
        self.nodes_lst.append( new_node )
        self.nodes_map[node_number] = new_node
        self.file_lineno += 1

      # numbering the fluent volume cells
      nnode = 1
      for nd in self.nodes_lst:
        nd.fl_number = nnode
        nnode += 1
      #
      rgx = re.compile( re_end_node )
      tokens = rgx.search( self.file_lines[ self.file_lineno ] )
      if tokens:
        self.file_lineno += 1
        return True
      else:
        return False
    else:
      return False

  ##----------------------------------------------------------------------------
  def __gmsh_read_cells( self ):
    rgx = re.compile( re_begin_elm )
    tokens = rgx.search( self.file_lines[ self.file_lineno ] )
    if tokens:
      self.file_lineno += 1
      fields = string.split( self.file_lines[ self.file_lineno ] )
      num_cells = int( fields[0] )
      self.file_lineno += 1

      # the allowed type as function of the dimension
      if self.mesh_dim == 2:
        self.cell_type[1] = []
        self.cell_type[2] = []
        self.cell_type[3] = []
      else:
        self.cell_type[2] = []
        self.cell_type[3] = []
        self.cell_type[4] = []
        self.cell_type[5] = []
        self.cell_type[6] = []
        self.cell_type[7] = []
      #
      for j in range( 0, num_cells ):
        fields = string.split( self.file_lines[ self.file_lineno ] )

        elm_number      = int( fields[0] )
        elm_type        = int( fields[1] )
        reg_phys        = int( fields[2] )
        reg_elem        = int( fields[3] )
        number_of_nodes = int( fields[4] )

        if number_of_nodes > 1:
          nodes_lst = []
          for i in range( 0, number_of_nodes ):
            nodes_lst.append( self.nodes_map[ int( fields[5+i] ) ] )

          if elm_type == 1 and self.mesh_dim == 3:
            self.file_lineno += 1
            continue

          new_cell = cell( elm_number, elm_type, reg_phys,
                            reg_elem, tuple( nodes_lst ) )

          self.cells_map[ elm_number ] = new_cell

          cell_type_lst = self.cell_type[ elm_type ]
          cell_type_lst.append( new_cell )

        self.file_lineno += 1

      # numbering the fluent volume cells
      ncell = 1
      if self.mesh_dim == 2:
        for cl in self.cell_type[2]:
          cl.fl_number = ncell
          ncell += 1
        for cl in self.cell_type[3]:
          cl.fl_number = ncell
          ncell += 1
      else:
        for cl in self.cell_type[4]:
          cl.fl_number = ncell
          ncell += 1
        for cl in self.cell_type[5]:
          cl.fl_number = ncell
          ncell += 1
        for cl in self.cell_type[6]:
          cl.fl_number = ncell
          ncell += 1
        for cl in self.cell_type[7]:
          cl.fl_number = ncell
          ncell += 1
      #
      rgx = re.compile( re_end_elm )
      tokens = rgx.search( self.file_lines[ self.file_lineno ] )
      if tokens:
        self.file_lineno += 1
        return True
      else:
        return False
    else:
      return False

  ##----------------------------------------------------------------------------
  def __gmsh_read_per( self ):
    rgx = re.compile( re_begin_per )
    tokens = rgx.search( self.file_lines[ self.file_lineno ] )
    if tokens:
      self.file_lineno += 1
      fields = string.split( self.file_lines[ self.file_lineno ] )
      num_per = int( fields[0] )
      self.file_lineno += 1
      for j in range( 0, num_per ):
        fields = string.split( self.file_lines[ self.file_lineno ] )
        id_number = int( fields[0] )
        master = int( fields[1] )
        slave  = int( fields[2] )
        self.periodic_map[ master ] = slave
        self.file_lineno += 1

      rgx = re.compile( re_end_per )
      tokens = rgx.search( self.file_lines[ self.file_lineno ] )
      if tokens:
        self.file_lineno += 1
        return True
      else:
        return False
    else:
      return False

  ##----------------------------------------------------------------------------
  def read_gmsh( self, filename ):
    fin=open( filename, 'r' )

    try:
      self.file_lines = fin.readlines()
    finally:
      fin.close()

    self.file_lineno = 0

    while self.file_lineno < len( self.file_lines ):

      line = self.file_lines[ self.file_lineno ]
      print "-> ", line,

      if len( re.sub(r'\s', '', line ) ) == 0: # empty line ?
        self.file_lineno += 1
        continue

      if self.__gmsh_read_nodes():
        continue
      if self.__gmsh_read_cells():
        continue
      if self.__gmsh_read_per():
        continue
      else:
        exit_error( "at line: " + line )

  ##----------------------------------------------------------------------------
  def write_gmsh( self, filename ):
    #
    # IMPORTANT NOTE: the numbering is unique for each element type
    # may be repeated for elements with different type
    #
    fout = open( filename, 'w' )

    fout.write( "$NOD\n" )
    fout.write( "%i\n" % len( self.nodes_map ) )
    for nd in self.nodes_map.values():
      fout.write( "%i\t%.8f\t%.8f\t%.8f\n" %
                  ( nd.node_number, nd.point[0], nd.point[1], nd.point[2] ) )
    fout.write( "$ENDNOD\n" )

    fout.write( "$ELM\n" )
    fout.write( "%i\n" % len( self.cells_map ) )

    for cl in self.cells_map.values():
      nn = len( cl.nodes_lst )
      fout.write( "%i\t%i\t%i\t%i\t%i" % ( cl.elm_number, cl.elm_type,
                                           cl.reg_phys, cl.reg_elem, nn ) )
      for nd in cl.nodes_lst:
        fout.write( "\t%i" % nd.node_number )
      fout.write( "\n" )

    fout.write( "$ENDELM\n" )

    # non-standard
    if len( self.periodic_map ) != 0:
      fout.write( "$PERNODES\n" )
      fout.write( "%i\n" % len( self.periodic_map ) )

      i = 1
      for ( k, v ) in self.periodic_map.iteritems():
        fout.write( "%i\t%i\t%i\n" % ( i+1, k, v ) )
        i += 1

      fout.write( "$ENDPERNODES\n" )

    fout.close()

  ##----------------------------------------------------------------------------
  def remove_lines( self ):
    for cell in self.cells_map.values():
      if cell.elm_type == 1:
        del self.cells_map[ cell.elm_number ]
    cell_id = 0
    for cell in self.cells_map.values():
      cell.elm_number = cell_id
      cell_id += 1

  ##----------------------------------------------------------------------------
  def write_fluent( self, exp_vol_id, export_srfs, filename, scale = 1.0 ):

    boundary_faces = {}
    faces_cnx_map = {}

    # Populate boundary_faces and faces_cnx_map
    for cell in self.cells_map.values():

      if self.mesh_dim == 3:

        if cell.elm_type == 2:
          boundary_faces[ cell.face_hash(0) ] = cell

        elif cell.elm_type == 3:

          boundary_faces[ cell.face_hash(0) ] = cell

        elif cell.elm_type >= 4 and cell.elm_type <= 7:
          nfc = len( cell.faces_hash )
          if cell.elm_type == 5: assert( nfc == 6 )
          for f in range( 0, cell.num_faces() ):
            if faces_cnx_map.has_key( cell.face_hash(f) ):
              fcnx = faces_cnx_map[ cell.face_hash(f) ]
              assert( fcnx.left_cell == None and fcnx.left_face == None )
              fcnx.left_cell = cell
              fcnx.left_face = f
            else:
              faces_cnx_map[ cell.face_hash(f) ] = face_cnx( cell, f )

        else:
          exit_error( "write_fluent: invalid cell type 3D" )

      else:

        if cell.elm_type == 1:
          boundary_faces[ cell.face_hash(0) ] = cell

        elif cell.elm_type == 2 or cell.elm_type == 3:

          for f in range( 0, cell.num_faces() ):
            if faces_cnx_map.has_key( cell.face_hash(f) ):
              fcnx = faces_cnx_map[ cell.face_hash(f) ]
              assert( fcnx.left_cell == None and fcnx.left_face == None )
              fcnx.left_cell = cell
              fcnx.left_face = f
            else:
              faces_cnx_map[ cell.face_hash(f) ] = face_cnx( cell, f )

        else:
          exit_error( "write_fluent: invalid cell type 2D" )

    # Update BC type
    for fcnx in faces_cnx_map.values():
      if fcnx.left_cell == None:
        hash  = fcnx.right_cell.face_hash( fcnx.right_face )
        fcnx.bc_type = boundary_faces[ hash ].reg_elem

    # CHECKING MAP and count boundaries
    INTERIOR_SURFACE = 0

    bc_type_cnt = {}
    bc_type_cnt[ INTERIOR_SURFACE ] = 0

    for srf in export_srfs:
      bc_type_cnt[ srf.gmsh_bc_type ] = 0

    for fcnx in faces_cnx_map.values():
      assert( ( fcnx.bc_type != None and fcnx.left_cell == None )
              or
              ( fcnx.bc_type == None and fcnx.left_cell != None ) )

      if fcnx.bc_type != None:
        bc_type_cnt[fcnx.bc_type] = bc_type_cnt[fcnx.bc_type] + 1
      else:
        bc_type_cnt[INTERIOR_SURFACE] = bc_type_cnt[INTERIOR_SURFACE] + 1

    num_nodes = len( self.nodes_map )

    #---------------------------------------------------------------------------
    file = open( filename, "w")

    file.write( "(0 \"Mesh converted by gmsh2gambit_nD\")\n" )
    file.write( "(0 \"Joao Henriques, IDMEC/IST, 2012\")\n\n" )
    file.write( "(0 \"DIMENSION:\")\n" )
    file.write( "(2 %i)\n\n" % self.mesh_dim )

    ##-------------------------------------------------------------
    file.write( "(0 \"NODES:\")\n" )
    file.write( "(10(0 1 %x 1 3))\n" % num_nodes  )
    file.write( "(10(1 1 %x 1 3)(\n" % num_nodes  )

    for node in self.nodes_lst:
      if self.mesh_dim == 3:
        #~ x = node.point[0] * scale
        #~ q = node.point[1] * scale
        #~ r = node.point[2] * scale
        #~ theta = q / r
        #~ y = r * math.sin( theta )
        #~ z = r * math.cos( theta )

        #~ file.write( "%.9f\t%.9f\t%.9f\n" % ( x, y, z ) )

        file.write( "%.9f\t%.9f\t%.9f\n" % ( node.point[0] * scale,  \
                                              node.point[1] * scale,  \
                                              node.point[2] * scale ) )
      else:
        file.write( "%.9f\t%.9f\n" % ( node.point[0] * scale,  \
                                        node.point[1] * scale ) )

    file.write( "))\n\n" )

    ##-------------------------------------------------------------
    file.write( "(0 \"NUMBER OF FACES:\")\n" )
    file.write( "(13 (0 1 %x 0))\n\n" % len( faces_cnx_map )  )

    ##-------------------------------------------------------------
    has_periodic_faces = False

    fi = 1 # first index
    li = 0 # last index

    for srf in export_srfs:
      if srf.fl_bc_type[0] == bc_type_periodic[0]:

        has_periodic_faces = True

        fi = li + 1
        li = fi + bc_type_cnt[ srf.gmsh_bc_type ] - 1

        per_fi = fi
        per_li = li

        file.write( "(0 \"PERIODIC FACES %s:\")\n" % srf.fl_name )
        file.write( "(13 (%x %x %x %x 0)(\n" % ( srf.fl_id, fi, li,
                                                  bc_type_periodic[0] ) )

        master_cnt2slave = {}
        periodic_cnt = 1

        for fcnx in faces_cnx_map.values():
          if fcnx.bc_type == srf.gmsh_bc_type:

            nd_lst = []
            nodes = fcnx.right_cell.face_nodes( fcnx.right_face )
            file.write( "%x" % len( nodes ) )
            for nd in nodes:
              file.write( " %x" %  nd.fl_number )
              nd_lst.append( self.periodic_map[ nd.node_number ] )
            file.write( " %x %x\n" % ( fcnx.right_cell.fl_number, 0 ) )

            hash = hash_func( tuple( nd_lst ) )
            master_cnt2slave[ periodic_cnt ] = hash
            periodic_cnt += 1

        file.write( "))\n\n" )

    ##-------------------------------------------------------------
    for srf in export_srfs:
      if srf.fl_bc_type[0] == bc_type_periodic_shadow[0]:

        fi = li + 1
        li = fi + bc_type_cnt[ srf.gmsh_bc_type ] - 1

        shadow_face_indx = {}
        ii = fi

        file.write( "(0 \"SHADOW FACES %s:\")\n" % srf.fl_name )
        file.write( "(13 (%x %x %x %x 0)(\n" % ( srf.fl_id, fi, li,
                                                bc_type_periodic_shadow[0] ) )

        slave_hash2cnt = {}

        for fcnx in faces_cnx_map.values():
          if fcnx.bc_type == srf.gmsh_bc_type:

            nodes = fcnx.right_cell.face_nodes( fcnx.right_face )
            file.write( "%x" % len( nodes ) )
            for nd in nodes:
              file.write( " %x" %  nd.fl_number )
            file.write( " %x %x\n" % ( fcnx.right_cell.fl_number, 0 ) )

            hash = fcnx.right_cell.face_hash( fcnx.right_face )
            slave_hash2cnt[ hash ] = periodic_cnt
            periodic_cnt += 1

        file.write( "))\n\n" )

    ##-------------------------------------------------------------
    if has_periodic_faces == True:
      file.write( "(0 \"PERIODIC FACES PAIRS:\")\n" )
      file.write( "(18 (%x %x 3 4)(\n" % ( per_fi, per_li ) )

      for ( mst_cnt, slv_hash ) in master_cnt2slave.iteritems():
        file.write( "%x %x\n" % ( mst_cnt, slave_hash2cnt[ slv_hash ] ) )
        per_fi += 1

      file.write( "))\n\n" )

    ##-------------------------------------------------------------
    interior_surfs = 0

    for srf in export_srfs:

      if srf.fl_bc_type[0] == bc_type_periodic[0] or \
         srf.fl_bc_type[0] == bc_type_periodic_shadow[0]:
        continue

      if srf.fl_bc_type[0] == bc_type_interior[0]:
        interior_surfs += 1
        continue

      fi = li + 1
      li = fi + bc_type_cnt[ srf.gmsh_bc_type ] - 1

      file.write( "(0 \"SURFACE %s:\")\n" % srf.fl_name )
      file.write( "(13 (%x %x %x %x 0)(\n" % ( srf.fl_id, fi, li,
                                                srf.fl_bc_type[0] ) )

      for fcnx in faces_cnx_map.values():
        if fcnx.bc_type == srf.gmsh_bc_type:

          nodes = fcnx.right_cell.face_nodes( fcnx.right_face )
          file.write( "%x" % len( nodes ) )
          for nd in nodes:
            file.write( " %x" %  nd.fl_number )
          file.write( " %x %x\n" % ( fcnx.right_cell.fl_number, 0 ) )

      file.write( "))\n\n" )

    ##-------------------------------------------------------------
    assert( interior_surfs == 1 ) # for now we only allow one interior

    for srf in export_srfs:

      if srf.fl_bc_type[0] != bc_type_interior[0]:
        continue

      fi = li + 1
      li = fi + bc_type_cnt[ srf.gmsh_bc_type ] - 1

      file.write( "(0 \"SURFACE %s:\")\n" % srf.fl_name )
      file.write( "(13 (%x %x %x %x 0)(\n" % ( srf.fl_id, fi, li,
                                                srf.fl_bc_type[0] ) )

      for fcnx in faces_cnx_map.values():
        if fcnx.bc_type == None:

          nodes = fcnx.right_cell.face_nodes( fcnx.right_face )
          file.write( "%x" % len( nodes ) )
          for nd in nodes:
            file.write( " %x" %  nd.fl_number )
          file.write( " %x %x\n" % ( fcnx.right_cell.fl_number,
                                      fcnx.left_cell.fl_number ) )

      file.write( "))\n\n" )

    ##-------------------------------------------------------------
    file.write( "(0 \"CELLS:\")\n" )


    if self.mesh_dim == 2:
      nv_2 = len( self.cell_type[2] )
      nv_3 = len( self.cell_type[3] )
      nv_g = nv_2 + nv_3

      file.write( "(12 (0 1 %x 0))\n" % nv_g )
      if nv_2 > 0 and nv_3 == 0:
        file.write( "(12 (%x %x %x 1 1))\n" % ( exp_vol_id.fl_id, 1, nv_2 ) )
      elif nv_3 > 0 and nv_2 == 0:
        file.write( "(12 (%x %x %x 1 3))\n" % ( exp_vol_id.fl_id, 1, nv_3 ) )
      else:
        file.write( "(12 (%x %x %x 1 0)(\n" % ( exp_vol_id.fl_id, 1, nv_g ) )
        col = 0

        for i in range( 0, nv_2 ):
          file.write( " 1" )
          col += 1
          if col % 9 == 0: file.write( "\n" )

        for i in range( 0, nv_3 ):
          file.write( " 3" )
          col += 1
          if col % 9 == 0: file.write( "\n" )

        if col % 9 != 0:
          file.write( "\n" )

        file.write( "))\n\n" )

    else:
      nv_4 = len( self.cell_type[4] )
      nv_5 = len( self.cell_type[5] )
      nv_6 = len( self.cell_type[6] )
      nv_7 = len( self.cell_type[7] )
      nv_g = nv_4 + nv_5 + nv_6 + nv_7

      file.write( "(12 (0 1 %x 0))\n" % nv_g )
      file.write( "(12 (%x %x %x 1 0)(\n" % ( exp_vol_id.fl_id, 1, nv_g ) )
      col = 0

      for i in range( 0, nv_4 ):
        file.write( " 2" )
        col += 1
        if col % 9 == 0: file.write( "\n" )

      for i in range( 0, nv_5 ):
        file.write( " 4" )
        col += 1
        if col % 9 == 0: file.write( "\n" )

      for i in range( 0, nv_6 ):
        file.write( " 6" )
        col += 1
        if col % 9 == 0: file.write( "\n" )

      for i in range( 0, nv_7 ):
        file.write( " 5" )
        col += 1
        if col % 9 == 0: file.write( "\n" )

      if col % 9 != 0:
        file.write( "\n" )

      file.write( "))\n\n" )

    ##-------------------------------------------------------------
    file.write( "(0 \"ZONES:\")\n" )

    # GAMBIT/FLUENT ERROR:
    # The id here is in decimal format but in hex in all other declarations
    file.write( "(45 (%i fluid %s)())\n" % ( exp_vol_id.fl_id,
                                              exp_vol_id.fl_name ) )
    for srf in export_srfs:
      file.write( "(45 (%i %s %s)())\n" % ( srf.fl_id, srf.fl_bc_type[1],
                                             srf.fl_name ) )

    file.close()

#==EOF==========================================================================