/***************************************************************************
                            mesh generation code
                            --------------------
    First version   [ 001 ] : Oct 17, 2001, 20:02
    Current version [ 088 ] : Jan 05, 2002, 14:15
    copyright               : (C) 2001 by Joao Carlos de Campos Henriques
    email                   : jcch@popsrv.ist.utl.pt
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t_mesh2d.h"

namespace mesh_2d
{
   const int FAC2D[3][2] = { { 1,2 },{ 2,0 },{ 0,1 } }; // face definition
   const int SUC2D[3][2] = { { 2,1 },{ 0,2 },{ 1,0 } }; // rotl
}

using namespace mesh_2d;

void mesh_2d::check_stream( FILE *stream, char *filename )
{
  if( stream == NULL )
  {
    printf( "Error opening file: '%s'", filename );
    exit(-1);
  }
}

mesh2d_base::mesh2d_base( node2d* (*ndalloc)(int), cell2d* (*clalloc)(int) ,
                          link2d* (*lkalloc)(int), edge2d* (*edalloc)() )
{
   node_alloc = ( ndalloc ? ndalloc : def_node_alloc );
   cell_alloc = ( clalloc ? clalloc : def_cell_alloc );
   link_alloc = ( lkalloc ? lkalloc : def_link_alloc );
   edge_alloc = ( edalloc ? edalloc : def_edge_alloc );
   cur_node_id = cur_cell_id = cur_link_id = 0;
   store_version = 0;
}

mesh2d_base::~mesh2d_base() {}

node2d* mesh2d_base::get_first_node()
{
   nodes_iterator = mesh_nodes.begin();
   return( nodes_iterator != mesh_nodes.end() ? *nodes_iterator : 0 );
}

node2d* mesh2d_base::get_next_node()
{
   nodes_iterator++;
   return( nodes_iterator != mesh_nodes.end() ? *nodes_iterator : 0 );
}

cell2d* mesh2d_base::get_first_cell()
{
   cells_iterator = mesh_cells.begin();
   return( cells_iterator != mesh_cells.end() ? *cells_iterator : 0 );
}

cell2d* mesh2d_base::get_next_cell()
{
   cells_iterator++;
   return( cells_iterator != mesh_cells.end() ? *cells_iterator : 0 );
}

node2d* mesh2d_base::def_node_alloc( int _id )
{
   //assert( false );
   return new node2d( _id );
}

cell2d* mesh2d_base::def_cell_alloc( int _id )
{
   return new cell2d( _id );
}

link2d* mesh2d_base::def_link_alloc( int _id )
{
   return new link2d( _id );
}

edge2d* mesh2d_base::def_edge_alloc()
{
   return new edge2d;
}

void
mesh2d_base::get_mesh_properties( int* nodes, int* cells )
{
   *nodes = mesh_nodes.size();
   *cells = mesh_cells.size();
}

void
mesh2d_base::load_node( FILE *stream, node2d** node_table, cell2d** /*cell_table*/ )
{
   int nid;

   tfread( nid, stream );
   node2d* nd = node_table[nid];

   tfread( nd->p.x, stream );
   tfread( nd->p.y, stream );
   tfread( nd->bc_type, stream );
   tfread( nd->bc_index, stream );
   tfread( nd->bc_surface, stream );
   nd->id = nid;
}

void
mesh2d_base::save_node( FILE *stream, node2d* nd )
{
  tfwrite( nd->id, stream );
  tfwrite( nd->p.x, stream );
  tfwrite( nd->p.y, stream );
  tfwrite( nd->bc_type, stream );
  tfwrite( nd->bc_index, stream );
  tfwrite( nd->bc_surface, stream );
}

void
mesh2d_base::load_cell( FILE *stream, node2d** node_table, cell2d** cell_table )
{
   int f, cid, nid;
   face2d* face;
   cell2d* cl;

   tfread( cid, stream );
   cl = cell_table[cid];
   cl->id = cid;

   tfread( cl->xc, stream );
   tfread( cl->yc, stream );
   tfread( cl->rc, stream );

   for( f = 0 ; f < 3 ; f++ )
   {
      tfread( nid, stream );
      face = &cl->face[f];
      face->node = node_table[nid];
   }

   cell_normals( cl );
}

void
mesh2d_base::save_cell( FILE *stream, cell2d* cl )
{
   int f, nid;
   face2d* face;

   tfwrite( cl->id, stream );
   tfwrite( cl->xc, stream );
   tfwrite( cl->yc, stream );
   tfwrite( cl->rc, stream );

   for( f = 0 ; f < 3 ; f++ )
   {
      face = &cl->face[f];
      nid  = face->node->id;
      tfwrite( nid, stream );
   }
}

void
mesh2d_base::load( FILE *stream, void(*progress)(int) )
{
   int num_nodes, num_cells;
   int i, total_work, work, cur = 0;

   tfread( store_version, stream );
   tfread( num_nodes, stream );
   tfread( num_cells, stream );
   total_work = num_cells + num_nodes;

   if( store_version > 1 )
      THROW__X( "mesh2d_base::load invalid version number.\n" );

   node2d** node_table = new node2d*[num_nodes+1];
   for( i = 1 ; i <= num_nodes ; i++ )
   {
      node_table[i] = node_alloc( ++cur_node_id );
      mesh_nodes.push_back( node_table[i] );
   }
   node_table[0] = 0;

   cell2d** cell_table = new cell2d*[num_cells+1];
   for( i = 1 ; i <= num_cells ; i++ )
   {
      cell_table[i] = cell_alloc( ++cur_cell_id );
      mesh_cells.insert( cell_table[i] );
   }
   cell_table[0] = 0;

   for( i = 0 ; i < num_nodes ; i++ )
   {
      load_node( stream, node_table, cell_table );
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
   }

   for( i = 0 ; i < num_cells ; i++ )
   {
      load_cell( stream, node_table, cell_table );
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
   }

   delete[] node_table;
   delete[] cell_table;

   recover_adjacent_links();
}

void
mesh2d_base::save( FILE *stream, void(*progress)(int) )
{
   int num_cells, num_nodes;
   int cur_id, total_work, work, cur = 0;

   node2d* nd;
   cell2d* cl;

   store_version = 1;

   get_mesh_properties( &num_nodes, &num_cells );
   total_work = num_cells + num_nodes;

   tfwrite( store_version, stream );
   tfwrite( num_nodes, stream );
   tfwrite( num_cells, stream );

   nd = get_first_node();
   cur_id = 0;
   while( nd )
   {
      nd->id = ++cur_id; // renumbering to load correctly
      save_node( stream, nd );
      nd = get_next_node();
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
   }

   cl = get_first_cell();
   cur_id = 0;
   while( cl )
   {
      cl->id = ++cur_id;
      save_cell( stream, cl );
      cl = get_next_cell();
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
   }
}

void
mesh2d_base::save_gmsh( FILE *stream, void(*progress)(int) )
{
   int num_cells, num_nodes, bc_faces = 0;
   int cur_id, total_work, work, cur = 0;
   face2d* face;

   node2d* nd;
   cell2d* cl;

   get_mesh_properties( &num_nodes, &num_cells );

  // count boundary faces
   cl = get_first_cell();
   while( cl )
   {
      for( int f_adj = 0 ; f_adj < 3 ; f_adj++ )
      {
        if( cl->face[ f_adj ].adj == 0 )
          bc_faces++;
      }
      cl = get_next_cell();
   }

   total_work = num_cells + bc_faces + num_nodes;

   fprintf( stream, "$NOD\n" );
   fprintf( stream, "%i\n", num_nodes );

   cur_id = 0;
   nd = get_first_node();
   while( nd )
   {
      nd->id = cur_id++; // renumbering to load correctly
      fprintf( stream, "%6i  % .12f  % .12f  % .1f\n", nd->id, nd->p.x, nd->p.y, 0.0 );
      nd = get_next_node();
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
   }

   fprintf( stream, "$ENDNOD\n" );
   fprintf( stream, "$ELM\n" );
   fprintf( stream, "%i\n", num_cells + bc_faces );

   cl = get_first_cell();
   cur_id = 0;
   while( cl )
   {
      face = &cl->face[2];
      int nid0 = face->node->id;
      face = &cl->face[1];
      int nid1 = face->node->id;
      face = &cl->face[0];
      int nid2 = face->node->id;

      fprintf( stream, "%6i  2  0  0  3  %6i  %6i  %6i\n", cur_id++, nid0, nid1, nid2 );

      cl = get_next_cell();
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
   }

   cl = get_first_cell();
   while( cl )
   {
      for( int f_adj = 0 ; f_adj < 3 ; f_adj++ )
      {
         if( cl->face[ f_adj ].adj == 0 )
         {
            node2d* nd_adj_1 = cl->face[ FAC2D[f_adj][0] ].node;
            node2d* nd_adj_0 = cl->face[ FAC2D[f_adj][1] ].node;
            int nid0 = nd_adj_0->id;
            int nid1 = nd_adj_1->id;
            int bc_type = nd_adj_0->bc_type & nd_adj_1->bc_type;
            fprintf( stream, "%6i  1  0  %i  2  %6i  %6i\n", cur_id++, bc_type, nid0, nid1 );
         }
      }
      if( progress )
      {
         work = (++cur) * 100 / total_work;
         progress( work );
      }
      cl = get_next_cell();
   }
   fprintf( stream, "$ENDELM\n\n" );
}

void
mesh2d_base::recover_adjacent_links()
{
   static int num = 1;

   node2d *nd, *nf[2], *nd_adj_0, *nd_adj_1;
   cell2d *cl, *cl_adj;
   face2d *ln;
   int j, f, f_adj;
   bool face_founded, af_0, af_1;

   /*
    * clean actual lists
    */
   nd = get_first_node();
   while( nd )
   {
      nd->head = 0;
      nd = get_next_node();
   }

   cl = get_first_cell();
   while( cl )
   {
      for( f = 0 ; f < 3 ; f++ )
      {
         cl->face[f].next = 0;
         cl->face[f].adj = 0;
      }
      cl = get_next_cell();
   }
   /*
    * link cells to nodes
    */
   cl = get_first_cell();
   while( cl )
   {
      add_cell_to_nodes( cl );
      cl = get_next_cell();
   }
   /*
    * for each cell face find the adjacent face using the list of cells stored at nodes
    */
   cl = get_first_cell();
   while( cl )
   {
      for( f = 0 ; f < 3 ; f++ )
      {
         /*
          * Se a face já esta ligada passar à seguinte
          */
         if( cl->face[f].adj == 0 )
         {
            face_founded = false;

            nf[0] = cl->face[ FAC2D[f][0] ].node;
            nf[1] = cl->face[ FAC2D[f][1] ].node;

            for( j = 0 ; j < 2 && !face_founded ; j++ )
            {
               ln = node2d_get_first_face( nf[j] );
               while( ln && !face_founded )
               {
                  cl_adj = ln->cell;
                  ln = node2d_get_next_face( ln );
                  if( cl != cl_adj )
                  {
                     /*
                      * temos que encontrar a face correcta para ser ligada
                      */
                     for( f_adj = 0 ; f_adj < 3 ; f_adj++ )
                     {
                        nd_adj_0 = cl_adj->face[ FAC2D[f_adj][0] ].node;
                        nd_adj_1 = cl_adj->face[ FAC2D[f_adj][1] ].node;
                        /*
                         * este algoritmo não tem em conta a ordem dos nós
                         */
                        af_0 = ( nf[0] == nd_adj_0 || nf[0] == nd_adj_1 );
                        af_1 = ( nf[1] == nd_adj_0 || nf[1] == nd_adj_1 );

                        if( af_0 && af_1 )
                        {
                           cl->face[f].adj = &cl_adj->face[f_adj];
                           cl_adj->face[f_adj].adj = &cl->face[f];
                           face_founded = true;
                           //printf( "faces linked == %i\n", num );
                           num++;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }
      cl = get_next_cell();
   }
   adjacent_linked = true;
}

