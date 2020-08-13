/***************************************************************************
                            mesh generation code
                            --------------------
    First version   [ 001 ] : Oct 17, 2001, 20:02
    Current version [ 088 ] : Jan 05, 2002, 14:15
    copyright               : (C) 2001 by Joao Carlos de Campos Henriques
    email                   : jcch@popsrv.ist.utl.pt
 ***************************************************************************/

#ifndef T_MESH2D_H
#define T_MESH2D_H

#include <functional>
#include <list>
#include <stdio.h>
#include <set>
#include <math.h>
#include "efread.h"
#include "common.h"
#include <assert.h>

using namespace std;

namespace mesh_2d {

extern const int FAC2D[3][2];
extern const int SUC2D[3][2];

struct face2d;
struct cell2d;
struct link2d;

struct p2d
{
  double    x;
  double    y;

  p2d( double _x, double _y ) { x=_x; y=_y; }
  p2d() { x = y = 0.0; }
};

struct node2d
{
   p2d       p;
   face2d   *head;          // First cell owning the node
   int       id;

   int       bc_type;       // boundary condition at node
   int       bc_index;      // periodic node index
   int       bc_surface;    // surface id - used for dump

   double    param;
   double    h;             // spacing function
   double    degree;

   node2d( int _id = -1 )    { id = _id; head = 0; degree = param = h = 0.0; bc_type = bc_index = bc_surface = 0; }
   node2d( int _id, p2d *q ) { id = _id; head = 0; p = *q; degree = param = h = 0.0; bc_type = bc_index = bc_surface = 0; }
};

struct link2d
{
   node2d     *node;
   link2d     *next;
   link2d     *prev;
   face2d     *adj;
   double      angle;
   int         id;

   link2d( int _id ) { id = _id; }
};

struct face2d
{
   cell2d  *cell;          // pointer to the cell (DO NOT CHANGE)
   int      id;            // face number (DO NOT CHANGE)
   node2d  *node;          // node opposite to this face
   node2d  *mid_node;      // mid node of this face
   face2d  *adj;           // adjacent cell to this face
   face2d  *next;          // next cell owning the node
   int      orientation;
   p2d      normal;
};

struct cell2d
{
   face2d  face[3];
   double  xc;
   double  yc;
   double  rc;
   double  area;
   double  dh_dx;         // spacing function gradient
   double  dh_dy;
   bool    bad;
   int     id;
   double  max_length;

   cell2d( int _id )
   {
      for( int i = 0 ; i < 3 ; i++ )
      {
         face[i].cell = this;
         face[i].id = i;
         face[i].node = 0;
         face[i].mid_node = 0;
         face[i].adj = 0;
         face[i].orientation = 0;
         face[i].next = 0;
      }
      max_length = xc = yc = rc = area = 0.0;
      bad = true;
      id = _id;
      dh_dx = dh_dy = 0.0;
   }
};

struct front2d
{
   link2d  *local_first;
   link2d  *local_prev;
};

struct edge2d
{
   node2d     *n0;
   node2d     *n1;
   node2d     *new_node;
   face2d     *adj;
   cell2d     *ins_cell;
   double    h, theta, ro;
};

struct smaller_angl
{
   bool operator () ( link2d* a, link2d* b ) const
   {
      if( a->angle == b->angle )
         return( b->id < a->id );                         // secundary key
      else
         return( a->angle < b->angle );                   // primary key
   }
};

struct smaller_addr
{
   bool operator () ( link2d* a, link2d* b ) const
   {
      if( a->node->id == b->node->id )
         return( a->next->node->id < b->next->node->id ); // secundary key
      else
         return( a->node->id < b->node->id );             // primary key
   }
};

// edge2d must be found by it's nodes
//
struct smaller_eaddr
{
   bool operator () ( edge2d* a, edge2d* b ) const
   {
      if( a->n0->id == b->n0->id )
         return( a->n1->id < b->n1->id );                 // secundary key
      else
         return( a->n0->id < b->n0->id );                 // primary key
   }
};

struct smaller_eht
{
   bool operator () ( edge2d* a, edge2d* b ) const
   {
      double hrel = ( a->h - b->h ) / ( a->h + b->h );
      if( fabs( hrel ) <= 1E-12 )
      {
#if 0
         double trel = ( a->theta - b->theta ) / ( a->theta + b->theta );
         if( fabs( trel ) <= 1E-12 )
            return( a->ro < b->ro );                      // thirdary key
         else
            return( a->theta < b->theta );                // secundary key
#else
         double rrel = ( a->ro - b->ro ) / ( a->ro + b->ro );
         if( fabs( rrel ) <= 1E-12 )
            return( a->theta < b->theta );                // thirdary key
         else
            return( a->ro < b->ro );                      // secundary key
#endif
      }
      else
         return( a->h < b->h );                           // primary key
   }
};

struct smaller_cell
{
   bool operator () ( cell2d* a, cell2d* b ) const
   {
      return( b->id < a->id );
   }
};

typedef list<node2d*>                  node2d_list;
typedef set<cell2d*, smaller_cell>     cell2d_set;

class tri6_xda_interpolator
{
public:
  virtual void interpolate( node2d*, node2d*, node2d* ) = 0;
  virtual ~tri6_xda_interpolator() = 0;
};

class mesh2d_base
{
protected:


   //pthread_mutex_t       access_mutex;
   int                   store_version;

   cell2d_set            mesh_cells;
   node2d_list           mesh_nodes;

   bool                  adjacent_linked;

   node2d_list::iterator nodes_iterator;
   cell2d_set::iterator  cells_iterator;

   int                   cur_node_id;
   int                   cur_cell_id;
   int                   cur_link_id;

   node2d* (*node_alloc)( int );
   cell2d* (*cell_alloc)( int );
   link2d* (*link_alloc)( int );
   edge2d* (*edge_alloc)();

   static node2d* def_node_alloc( int );
   static cell2d* def_cell_alloc( int );
   static link2d* def_link_alloc( int );
   static edge2d* def_edge_alloc();

   void load_node( FILE*, node2d**, cell2d** );
   void save_node( FILE*, node2d* );

   void load_cell( FILE*, node2d**, cell2d** );
   void save_cell( FILE*, cell2d* );

   void recover_adjacent_links();
public:
   void lock_iterators();
   void unlock_iterators();

   node2d* get_first_node();
   node2d* get_next_node();

   cell2d* get_first_cell();
   cell2d* get_next_cell();

   void get_mesh_properties( int *nodes, int *cells );

   void load( FILE*, void(*progress)(int) = 0 );
   void save( FILE*, void(*progress)(int) = 0 );
   void save_gmsh( FILE*, void(*progress)(int) = 0 );

   mesh2d_base( node2d* (*ndalloc)(int) = 0, cell2d* (*clalloc)(int) = 0,
                link2d* (*lkalloc)(int) = 0, edge2d* (*edalloc)() = 0 );
   virtual ~mesh2d_base();
};

class fist2d : public mesh2d_base
{
protected:
   typedef set<link2d*, smaller_angl>     link2d_set_angl;
   typedef set<link2d*, smaller_addr>     link2d_set_addr;
   typedef set<edge2d*, smaller_eaddr>    edge2d_set_addr;

   link2d_set_addr   fist_front;          // front   nodes for the FIST triangulation
   link2d_set_angl   convex;              // convex  nodes for the FIST triangulation
   link2d_set_angl   reflex;              // reflex  nodes for the FIST triangulation
   link2d_set_addr   waiting;             // waiting nodes for the FIST triangulation
                                          // - reflex nodes have invalidate the
                                          //   triangulation of these nodes

   front2d           front;
   cell2d_set        bad_cells;           // stores the bad mesh cells
   node2d_list       mid_nodes;

   void     initialization();
   void     classify_store_link( link2d *cl );
   bool     reflex_nodes_in_cell( link2d *vim1, link2d *vi, link2d *vip1 );
   //bool     node_in_cone( link2d *vj, link2d *vim1, link2d *vi, link2d *vip1 );
   void     create_cell( link2d *vi, link2d *vip1, link2d *vim1 );
   void     erase_link_from_sets( link2d *v );
   void     attach_faces( face2d *fa, face2d *fb );

   void     add_length( link2d *l );
   bool     join_fronts();

   virtual  void  attach_bad_edges( cell2d * ) {}

public:
   void     begin_front();

   void     add_to_front( node2d* );
   void     add_to_front( double x, double y, double param, int bc_type, int bc_index = 0, int bc_surface = 0 );

   void     end_front();

   bool     fist_generation ();

            fist2d( node2d* (*ndalloc)(int) = 0, cell2d* (*clalloc)(int) = 0,
                    link2d* (*lkalloc)(int) = 0, edge2d* (*edalloc)() = 0 );
   virtual ~fist2d() {}
};

class mesh2d : public fist2d
{
protected:
   typedef set<edge2d*, smaller_eht>      edge2d_set_sort;

   cell2d_set      back_mesh;           // stores the background mesh cells

   double          x_centroid;
   double          y_centroid;

   edge2d_set_addr   frontal_edges;       // global front sorted by edge addr
                                          // - used to find an edge in global_front_size

   cell2d           *back_start;

   const char* dump_dir;
   bool dumping_enable;

   bool     check_dump_dir();

   double   interpolate( cell2d *bc, node2d *n2 );
   void     new_point( edge2d* be );
   bool     is_implicit( cell2d *c, int type );
   bool     is_frontal_face( face2d *f );
   bool     is_close_to_existing_node( cell2d *c, node2d *n );
   bool     is_point_allowed( cell2d *c, node2d *n );
   bool     is_close_to_frontal_face( cell2d *c, node2d *n );

   cell2d*  find_cell( cell2d_set* set, cell2d *c, node2d *n );
   cell2d*  find_back_cell( node2d *n );
   cell2d*  find_bad_cell( cell2d *c, node2d *n );

   bool     green_sibson( face2d *f1 );
   void     create_fist_front( cell2d *c, node2d *n );
   void     create_new_cell( edge2d *be, node2d *n );
   void     create_frontal_edges();
   void     detach_bad_edges( cell2d *c );
   virtual  void  attach_bad_edges( cell2d *c );

   void     centroide();
   void     create_back_mesh();
   void     implicit_cells();
   bool     make_delaunay( cell2d_set* );

public:
   bool     mesh_generation ();
   bool     is_delaunay();

   void     set_dump_dir( const char* );
   void     set_dumping( bool b ) { dumping_enable = b; }
   void     dump_bad_cells();
   void     dump_mesh();
   void     dump_links();
   void     dump_back();
   void     dump_frontal_edges();
   void     dump_cell( cell2d*, cell2d*, node2d *n );
   void     testing_mesh();
   void     smooth();

            mesh2d( node2d* (*ndalloc)(int) = 0, cell2d* (*clalloc)(int) = 0,
                    link2d* (*lkalloc)(int) = 0, edge2d* (*edalloc)() = 0);
   virtual ~mesh2d() {}

   cell2d_set&  get_mesh_cells_set() { return mesh_cells; }
   node2d_list& get_mesh_nodes_list() { return mesh_nodes; }
   node2d_list& get_mesh_mid_nodes_list() { return mid_nodes; }

   void convert_to_tri6( tri6_xda_interpolator& );
   void save_tri6_xda( FILE *stream );
};

inline void
add_cell_to_nodes( cell2d *cl )
{
   face2d *ln;
   node2d *nd;
   int     fc;
   for( fc = 0 ; fc < 3 ; fc++ )
   {
      ln = &cl->face[fc];
      nd = ln->node;
      ln->next = nd->head;
      nd->head = ln;
   }
}

inline void
remove_cell_from_nodes( cell2d *cl )
{
   face2d *ln_1, *ln_2;
   node2d *nd_1;
   int   fc;
   for( fc = 0 ; fc < 3 ; fc++ )
   {
      ln_1 = &cl->face[fc];
      nd_1 = ln_1->node;
      if( nd_1->head == ln_1 )
      {
         nd_1->head = ln_1->next;
         ln_1->next = 0;
      }
      else
      {
         ln_2 = ln_1->node->head;
         while( ln_2->next )
         {
            if( ln_2->next == ln_1 )
            {
               ln_2->next = ln_1->next;
               ln_1->next = 0;
               break;
            }
            else
            {
               ln_2 = ln_2->next;
               if( ln_2->next == 0 )
                  THROW__X( "remove_cell_from_nodes -> cell not found." );
            }
         }
      }
   }
}

inline face2d*
node2d_get_first_face( node2d *nd )
{
   return nd->head;
}

inline face2d*
node2d_get_next_face( face2d *ln )
{
   return ln->next;
}


extern const double mesh_eps;
extern const double S2min_S2max;
extern const double S2max_S2min;
extern const double dist_factor;
extern const int succ[3];
extern const int pred[3];
void check_stream( FILE *stream, char *filename );

#include "t_mesh2d_aux_funcs.h"

}; // namespace __mesh2d
#endif

//***EOF************************************************************************
