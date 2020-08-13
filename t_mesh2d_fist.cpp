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

using namespace mesh_2d;

fist2d::fist2d( node2d* (*ndalloc)(int), cell2d* (*clalloc)(int),
                link2d* (*lkalloc)(int), edge2d* (*edalloc)() )
              : mesh2d_base( ndalloc, clalloc, lkalloc, edalloc )
{
   if( ( sizeof( link2d* ) - sizeof( unsigned long ) ) != 0 )
      THROW__X( "fist2d::fist2d: program assumes sizeof( pointer )"
                " == sizeof( int ), which is false..." );
   begin_front();
}

void
fist2d::add_length( link2d *l )
{
   double h = 0.5 * norm( l->node->p.x - l->next->node->p.x,
                          l->node->p.y - l->next->node->p.y );
   l->node->h += h;
   l->next->node->h += h;
}

void
fist2d::begin_front()
{
   front.local_first = front.local_prev = 0;
}

void
fist2d::add_to_front( double x, double y, double param, int bc_type, int bc_index, int bc_surface )
{
   node2d *n = node_alloc(0);
   n->p.x = x;
   n->p.y = y;

   n->param      = param;
   n->bc_type    = bc_type;
   n->bc_index   = bc_index;
   n->bc_surface = bc_surface;

   add_to_front( n );
}

void
fist2d::add_to_front( node2d *n )
{
   n->id = ++cur_node_id;
   mesh_nodes.push_back( n );

   link2d *l = link_alloc( ++cur_link_id );
   l->node = n;
   l->prev = front.local_prev;
   l->next = 0;
   l->adj = 0;

   if( front.local_prev != 0 )
   {
      front.local_prev->next = l;
      fist_front.insert( front.local_prev );
      add_length( front.local_prev );
   }
   if( front.local_first == 0 )
      front.local_first = l;
   front.local_prev = l;
}

void
fist2d::end_front()
{
   front.local_first->prev = front.local_prev;
   front.local_prev->next = front.local_first;
   fist_front.insert( front.local_prev );
   add_length( front.local_prev );
}

/***********************************************************************
   FIST code
 ***********************************************************************/

void
fist2d::classify_store_link( link2d *cl )
{
   double ax, ay, bx, by, x, y;

   ax = cl->prev->node->p.x - cl->node->p.x;
   ay = cl->prev->node->p.y - cl->node->p.y;

   bx = cl->next->node->p.x - cl->node->p.x;
   by = cl->next->node->p.y - cl->node->p.y;

   y = ax * by - ay * bx;
   x = ax * bx + ay * by;

   cl->angle = atan2pi( y, x );

   if( cl->angle < M_PI )
      convex.insert( cl );
   else
      reflex.insert( cl );
}

void
fist2d::initialization()
{
   link2d_set_addr::iterator it = fist_front.begin();
   link2d *l;

   while( it != fist_front.end() )
   {
      l = *it;
      classify_store_link( l );
      it++;
   }
   fist_front.clear();
}

bool
fist2d::reflex_nodes_in_cell( link2d *vim1, link2d *vi, link2d *vip1 )
{
   link2d_set_angl::iterator it = reflex.begin();
   link2d* cl = *it;
   int a, b, c;

   while( it != reflex.end() )
   {
      cl = *it;

      if( vi->node != cl->node && vim1->node != cl->node && vip1->node != cl->node )
      {
         a = det( vi->node, vip1->node, cl->node );
         b = det( vip1->node, vim1->node, cl->node );
         c = det( vim1->node, vi->node, cl->node );

         if( a <= 0 && b <= 0 && c <= 0 )
            return true;
      }
      it++;
   }
   return false;
}

/*
bool
fist2d::node_in_cone( link2d *vj, link2d *vim1, link2d *vi, link2d *vip1 )
{
   int a, b;

   a = det( vi->node, vip1->node, vj->node );
   b = det( vi->node, vim1->node, vj->node );

   if( vi->angle <= M_PI )
      return( a < 0 && b > 0 );
   else
      return( a < 0 || b > 0 );
}
*/

void
fist2d::erase_link_from_sets( link2d *v )
{
   link2d_set_angl::iterator it;

   it = convex.find( v );
   if ( it != convex.end() )
   {
      convex.erase( it );
      return;
   }

   it = reflex.find( v );
   if ( it != reflex.end() )
   {
      reflex.erase( it );
      return;
   }

   it = waiting.find( v );
   if ( it != waiting.end() )
   {
      waiting.erase( it );
      return;
   }

   THROW__X( "fist2d::erase_link_from_sets: link2d not found in auxiliary sets." );
}

void
fist2d::attach_faces( face2d *fa, face2d *fb )
{
   if( fa != 0 ) fa->adj = fb;
   if( fb != 0 ) fb->adj = fa;
}

void
fist2d::create_cell( link2d *vi, link2d *vip1, link2d *vim1 )
{
   cell2d *c = cell_alloc( ++cur_cell_id );

   c->face[0].node = vi->node;
   c->face[1].node = vip1->node;
   c->face[2].node = vim1->node;

   circun_circle( c );
   bad_cells.insert( c );

   attach_bad_edges( c );  // used only in internal mesh generation

   erase_link_from_sets( vim1 );    // the link angles have changed
   erase_link_from_sets( vip1 );    // ..

   attach_faces( c->face+1, vim1->adj );
   attach_faces( c->face+2, vi->adj );

   if( vip1->next == vim1 )                    // three edge front
   {
      attach_faces( c->face+0, vip1->adj );    // the final link
      return;
   }

   vip1->prev = vim1;
   vim1->next = vip1;

   classify_store_link( vip1 );
   classify_store_link( vim1 );

   vim1->adj = c->face+0;
}

bool
fist2d::join_fronts()
{
   link2d_set_angl neighbours; // all candidates to the new triangulation

   link2d *vi   = *waiting.begin();
   link2d *vip1 = vi->next;
   link2d *vim1 = vi->prev;

   link2d *l0, *cl, *clm1;
   cell2d *cn;
   node2d *ni, *nip1, *nim1, *ncl, *nl0, *nl1, *nl2, mid( 0 );
   //
   // find all reflex nodes that are inside this "triangle" (vi)
   //
   link2d_set_angl::iterator itv, it = reflex.begin();
   link2d* lt;
   int a, b, c;

   ni   = vi->node;
   nip1 = vip1->node;
   nim1 = vim1->node;

   mid.p.x = 0.5 * ( ni->p.x + nip1->p.x );
   mid.p.y = 0.5 * ( ni->p.y + nip1->p.y );

   while( it != reflex.end() )
   {
      lt  = *it;
      nl1 = lt->node;
      if( ni != nl1 && nip1 != nl1 && nim1 != nl1 )
      {
         a = det( ni,   nip1, nl1 );
         b = det( nip1, nim1, nl1 );
         c = det( nim1, ni,   nl1 );

         if( a <= 0 && b <= 0 && c <= 0 )
            neighbours.insert( lt );
      }
      it++;
   }
   if( neighbours.empty() )
      THROW__X( "fist2d::join_fronts: waiting link without reflex inside." );
   //
   // Now that we have all neighbours, find one that creates a valid triangle.
   // For that, we check if another neighbour is inside the current triangle.
   //
   it = neighbours.begin();
   while( it != neighbours.end() )
   {
      cl  = *it;   // new candidate to form a valid triangle
      ncl = cl->node;

      itv = neighbours.begin();
      while( itv != neighbours.end() )
      {
         lt  = *itv;
         nl0 = lt->prev->node;
         nl1 = lt->node;
         nl2 = lt->next->node;

         if( ni != nl1 && nip1 != nl1 && ncl != nl1 )
         {
            //
            // test point inside
            //
            a = det( ni,   nip1, nl1 );
            b = det( nip1, ncl,  nl1 );
            c = det( ncl,  ni,   nl1 );
            if( a <= 0 && b <= 0 && c <= 0 )
               break;   // found an invalid triangle
            //
            // test if both edges of lt intersect mid line
            // this is a case that only hapens on front merging
            //
            if( ncl != nl2 && segs_intersection( &mid, ncl, nl1, nl2 ) )
               break;   // found an invalid triangle

            if( ncl != nl0 && segs_intersection( &mid, ncl, nl0, nl1 ) )
               break;   // found an invalid triangle
         }
         itv++;
      }
      //
      // If search terminated sucessfully then create the new cell
      //
      if( itv == neighbours.end() )
      {
         // CHECK the possible error - a node with several link2d ...
         clm1 = cl->prev;
         cn = cell_alloc( ++cur_cell_id );

         cn->face[0].node = ni;
         cn->face[1].node = nip1;
         cn->face[2].node = ncl;

         circun_circle( cn );
         bad_cells.insert( cn );
         attach_bad_edges( cn );

         erase_link_from_sets( vi );      // the link angles have changed
         erase_link_from_sets( vip1 );    // ..
         erase_link_from_sets( cl );      // ..

         attach_faces( cn->face+2, vi->adj );

         l0 = link_alloc( ++cur_link_id );
         l0->node = cl->node;
         l0->adj  = cn->face+0;
         l0->prev = clm1;
         l0->next = vip1;

         clm1->next = l0;
         vip1->prev = l0;

         cl->prev = vi;
         vi->next = cl;
         vi->adj  = cn->face+1;

         classify_store_link( vi );
         classify_store_link( cl );
         classify_store_link( l0 );
         classify_store_link( vip1 );

         return true;
      }
      it++;
   }
   printf( "fist2d::join_fronts: failed to join fronts\n" );
   return false;
}

bool
fist2d::fist_generation()
{
   link2d *vi, *vim1, *vip1, *vim2, *vip2;
   bool rflx_in_cell;
   //bool vim1_in_cone, vip1_in_cone;

   if( fist_front.empty() )
   {
      printf( "fist2d::fist_generation: nothing to do! Empty front list.\n" );
      return false;
   }

   initialization();

   while( !convex.empty() )
   {
      while( !convex.empty() )
      {
         vi = *convex.begin();
         convex.erase( vi );

         vip1 = vi->next;
         vip2 = vip1->next;

         vim1 = vi->prev;
         vim2 = vim1->prev;

         // (1) vi is taken from the convex list
         // (2) vi, vip1, vim1 does not contain any reflex node
         // (3) vim1 in C( vi, vip1, vip2 ) and vip1 in C( vim2, vim1, vi )

         rflx_in_cell = reflex_nodes_in_cell( vim1, vi, vip1 );
         if( vip2 == vim1 && !rflx_in_cell )
            create_cell( vi, vip1, vim1 );
         else
         {
            //vim1_in_cone = node_in_cone( vim1, vi, vip1, vip2 );
            //vip1_in_cone = node_in_cone( vip1, vim2, vim1, vi );

            //if( !rflx_in_cell && ( !vim1_in_cone || !vip1_in_cone ) )
            //   printf("##################################################################################\n" );

            if( !rflx_in_cell ) //&& vim1_in_cone && vip1_in_cone )
               create_cell( vi, vip1, vim1 );
            else
               waiting.insert( vi );    // waiting list for angle change during mesh generation
         }
      }

      if( waiting.empty() )
         return false;
      else
      {
         // try to create a triangle with another front
         if( !join_fronts() )
            THROW__X( "fist2d::fist_generation: failed to join fronts." );
      }
   }
   return true;
}


//***EOF************************************************************************
