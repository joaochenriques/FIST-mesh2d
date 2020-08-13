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
   const double golden_ratio = 0.5 * ( sqrt( 5.0 ) - 1.0 );
   const double mesh_eps  = 1E-15;

   const double S2min_S2max = 1.0 - golden_ratio;
   const double S2max_S2min = 1.0 / S2min_S2max;
   const double impl_factor = 1.0 / golden_ratio;
   const double dist_factor = golden_ratio;

   const int succ[3] = { 1, 2, 0 };
   const int pred[3] = { 2, 0, 1 };

   const char mov_up[] = "\x1B[A";
}

using namespace mesh_2d;

/***********************************************************************
   Internal mesh generation code
 ***********************************************************************/

mesh2d::mesh2d( node2d* (*ndalloc)(int), cell2d* (*clalloc)(int),
                link2d* (*lkalloc)(int), edge2d* (*edalloc)() )
              : fist2d( ndalloc, clalloc, lkalloc, edalloc )
{
  back_start = 0;
  dump_dir = 0;
  dumping_enable = true;
}

void
mesh2d::set_dump_dir( const char *dir )
{
  dump_dir = dir;
}

bool
mesh2d::is_frontal_face( face2d *f )
{
   if( f->adj == 0 )
      return true;
   else
      return( !f->adj->cell->bad );
}

bool
mesh2d::is_close_to_frontal_face( cell2d *c, node2d *n )
{
   face2d *f;
   node2d *n0, *n1;
   double xm, ym, rm, sm;
   int i;

   for( i = 0 ; i < 3 ; i++ )
   {
      f = c->face+i;
      if( is_frontal_face( f ) )
      {
         n0 = succ_node( f );
         n1 = pred_node( f );

         xm = 0.5 * ( n0->p.x + n1->p.x );
         ym = 0.5 * ( n0->p.y + n1->p.y );
         rm = norm_sqr( n0->p.x - n1->p.x, n0->p.y - n1->p.y );
         sm = norm_sqr( n->p.x - xm, n->p.y - ym );

         if( sm < rm )
            return true;
      }
   }
   return false;
}

void
mesh2d::centroide()
{
   cell2d_set::iterator it;
   cell2d *cp;
   double sc, c2;

   x_centroid = y_centroid = sc = 0.0;
   for( it = bad_cells.begin() ; it != bad_cells.end() ; it++ )
   {
      cp = *it;
      c2 = sqr( cp->rc );
      x_centroid += cp->xc * c2;
      y_centroid += cp->yc * c2;
      sc += c2;
   }
   x_centroid /= sc;
   y_centroid /= sc;
}

void
mesh2d::create_back_mesh()
{
   cell2d_set::iterator itc, itb, ita;
   cell2d *cc, *cb, *ca;
   face2d *f;
   node2d *nn, *np, *ns;
   double nx, ny, meas;
   int i, fid;

   make_delaunay( &bad_cells );
   //
   // create back cells
   //
   for( itc = bad_cells.begin() ; itc != bad_cells.end() ; itc++ )
   {
      cc = *itc;
      cb = cell_alloc( ++cur_cell_id );
      cb->id = cc->id;
      for( i = 0 ; i < 3 ; i++ )
         cb->face[i].node = cc->face[i].node;
      back_mesh.insert( cb );
   }
   //
   // attach back cells faces
   //
   for( itb = back_mesh.begin(), itc = bad_cells.begin() ; itb != back_mesh.end() ; itb++, itc++ )
   {
      cb = *itb;
      cc = *itc;

      if( cb->id != cc->id )
         THROW__X( "mesh2d::create_back_mesh: back != cell mesh.\n" );

      for( i = 0 ; i < 3 ; i++ )
      {
         if( cb->face[i].adj != 0 )    // already attached
            continue;
         f = cc->face[i].adj;
         if( f != 0 )
         {
            fid = f->id;
            ita = back_mesh.find( f->cell );
            if( ita == back_mesh.end() )
               THROW__X( "mesh2d::create_back_mesh: back mesh cell not found.\n" );
            ca = *ita;
            attach_faces( cb->face+i, ca->face+fid );
         }
         else
            cb->face[i].adj = 0;

      }
   }
   //
   // evaluate cell gradient
   //
   for( itb = back_mesh.begin() ; itb != back_mesh.end() ; itb++ )
   {
      cb = *itb;

      for( i = 0 ; i < 3 ; i++ )
      {
         f = cb->face+i;

         nn = f->node;
         ns = succ_node( f );
         np = pred_node( f );

         nx = + np->p.y - ns->p.y;
         ny = - np->p.x + ns->p.x;

         cb->dh_dx += nn->h * nx;
         cb->dh_dy += nn->h * ny;
      }

      meas = cell_area( cb );
      cb->dh_dx /= meas;
      cb->dh_dy /= meas;
   }
}

double
mesh2d::interpolate( cell2d *bc, node2d *n2 )
{
   double X1 = bc->face[0].node->p.x;
   double Y1 = bc->face[0].node->p.y;
   double h1 = bc->face[0].node->h;

   double X2 = bc->face[2].node->p.x;
   double Y2 = bc->face[2].node->p.y;
   double h2 = bc->face[2].node->h;

   double X3 = bc->face[1].node->p.x;
   double Y3 = bc->face[1].node->p.y;
   double h3 = bc->face[1].node->h;

   double ar = ( X2*Y3 + X1*Y2 + X3*Y1 - Y1*X2 - Y2*X3 - Y3*X1 );
   double a1 = ( h1*(X2*Y3-X3*Y2) + h2*(X3*Y1-X1*Y3) + h3*(X1*Y2-X2*Y1) ) / ar;
   double a2 = ( h1*(Y2-Y3) + h2*(Y3-Y1) + h3*(Y1-Y2) ) / ar;
   double a3 = ( h1*(X3-X2) + h2*(X1-X3) + h3*(X2-X1) ) / ar;

   return( a1 + a2 * n2->p.x + a3 * n2->p.y );
}

void
mesh2d::new_point( edge2d *be )
{
   static const double sin60 = 0.5 * sqrt( 3.0 );

   node2d *p0 = be->n0;
   node2d *p1 = be->n1;
   node2d *n2 = be->new_node;

   double nx = p1->p.y - p0->p.y;
   double ny = p0->p.x - p1->p.x;
   double nn = norm( nx, ny );

   nx /= nn;
   ny /= nn;

   n2->p.x = 0.5 * ( p0->p.x + p1->p.x );
   n2->p.y = 0.5 * ( p0->p.y + p1->p.y );

   n2->p.x += 1E-6 * nn * nx;
   n2->p.y += 1E-6 * nn * ny;

   cell2d *cb = find_back_cell( n2 );

   double hm = interpolate( cb, n2 );
   double lm = hm / ( sin60 - 0.5 * ( cb->dh_dx * nx + cb->dh_dy * ny ) );

   double lmin = nn * sqrt( S2min_S2max - 0.25 );
   double lmax = nn * sqrt( S2max_S2min - 0.25 );

   lm = max( lmin, min( lm, lmax ) );
   n2->p.x += lm * nx;
   n2->p.y += lm * ny;
   n2->h  = hm;

   be->ins_cell = find_bad_cell( be->adj->cell, n2 );
}

cell2d*
mesh2d::find_cell( cell2d_set* st, cell2d *cl, node2d *n )
{
   // cl is the starting cell
   int i, a[3];

   while( cl != 0 )
   {
      a[0] = det( cl->face[1].node, cl->face[2].node, n );
      a[1] = det( cl->face[2].node, cl->face[0].node, n );
      a[2] = det( cl->face[0].node, cl->face[1].node, n );

      if( a[0] <= 0 && a[1] <= 0 && a[2] <= 0 )
         return cl;

      for( i = 0 ; i < 3 ; i++ )
         if( a[i] > 0 && cl->face[i].adj )
         {
            cl = cl->face[i].adj->cell;
            continue;
         }
      //
      // co-linear case, it is our last chance
      //
      for( i = 0 ; i < 3 ; i++ )
         if( a[i] == 0 && cl->face[i].adj )
         {
            cl = cl->face[i].adj->cell;
            continue;
         }
      break;
   }
   //
   // because the domain is not always convex we may need a global search
   //
   cell2d_set::iterator it = st->begin();
   while( it != st->end() )
   {
      cl = *it;

      a[0] = det( cl->face[1].node, cl->face[2].node, n );
      a[1] = det( cl->face[2].node, cl->face[0].node, n );
      a[2] = det( cl->face[0].node, cl->face[1].node, n );

      if( a[0] <= 0 && a[1] <= 0 && a[2] <= 0 )
         return cl;
      it++;
   }
   return 0;
}

cell2d*
mesh2d::find_back_cell( node2d *n )
{
   if( back_start == 0 )
      back_start = *back_mesh.begin();
   back_start = find_cell( &back_mesh, back_start, n );
   return back_start;
}

cell2d*
mesh2d::find_bad_cell( cell2d *cl, node2d *n )
{
   return find_cell( &bad_cells, cl, n );
}

void
mesh2d::detach_bad_edges( cell2d *c )
{
   edge2d_set_addr::iterator it;
   edge2d e, *a;
   int i;

   for( i = 0 ; i < 3 ; i++ )
   {
      e.n0 = succ_node( c->face+i );
      e.n1 = pred_node( c->face+i );

      it = frontal_edges.find( &e );
      if( it != frontal_edges.end() )
      {
         a = *it;
         a->adj = 0;
      }
   }
}

inline bool
inside_circuncircle( cell2d *c, node2d *n )
{
   double r = norm( c->xc - n->p.x, c->yc - n->p.y ) * ( 1.0 - 1E-8 );
   return( r <= c->rc );
}

void
attach_links( link2d *prv, link2d *nxt )
{
   prv->next = nxt;
   nxt->prev = prv;
}

void
mesh2d::create_fist_front( cell2d *c, node2d *n )
{
   //
   // it is assumed that c contains n
   //
   link2d_set_addr::iterator itl;
   link2d_set_addr untested;

   cell2d_set::iterator itt;
   cell2d_set trash;

   link2d *l[3], *cur, *l0, *l1, *prv, *nxt;
   face2d *adj0, *adj1, *face;
   node2d *nod0, *nod1;
   int i;

   if( fist_front.size() != 0 )
      THROW__X( "mesh2d::create_fist_front: fist_front.size() != 0\n" );

   for( i = 0 ; i < 3 ; i++ )
   {
      l[i] = link_alloc( ++cur_link_id );
      l[i]->node  = c->face[i].node;
   }
   for( i = 0 ; i < 3 ; i++ )
   {
      l[i]->next = l[ succ[i] ];
      l[i]->prev = l[ pred[i] ];
      l[i]->adj  = c->face[ pred[i] ].adj;
      untested.insert( l[i] );
   }
   trash.insert( c );

   while( !untested.empty() )
   {
      itl = untested.begin();
      cur = *itl;

      untested.erase( itl );
      face = cur->adj;

      if( face != 0 && inside_circuncircle( face->cell, n ) )
      {
         prv = cur->prev;
         nxt = cur->next;

         nod0 = cur->node;
         nod1 = face->node;

         adj0 = succ_face( face )->adj;
         adj1 = pred_face( face )->adj;

         itt = trash.find( face->cell );
         if( itt == trash.end() )
            trash.insert( face->cell );

         l0 = cur;
         l0->node = nod0;
         l0->adj = adj0;

         l1 = link_alloc( ++cur_link_id );
         l1->node = nod1;
         l1->adj = adj1;

         attach_links( prv, l0 );
         attach_links( l0, l1 );
         attach_links( l1, nxt );

         untested.insert( l0 );
         untested.insert( l1 );
         //
         // test unusual cases
         //
         if( prv->node == l0->next->node )
         {
            untested.erase( prv );
            untested.erase( l0 );
            attach_links( prv->prev, l0->next );
            delete prv;
            delete l0;
         }
         if( l1->node == nxt->next->node )
         {
            untested.erase( l1 );
            untested.erase( nxt );
            attach_links( l1->prev, nxt->next );
            delete l1;
            delete nxt;
         }
      }
      else
         fist_front.insert( cur );
   }
   //
   // delete all intersected cells
   //
   for( itt = trash.begin() ; itt != trash.end() ; itt++ )
   {
      c = *itt;
      for( i = 0 ; i < 3 ; i++ )
      {
         face = c->face[i].adj;
         if( face != 0 )
            face->adj = 0;
      }
      bad_cells.erase( c );
      detach_bad_edges( c );
      delete c;
   }
}

void
mesh2d::create_new_cell( edge2d *be, node2d *n )
{
   link2d_set_addr::iterator it;
   link2d *lf, *l0, *vim1, *vip1;

   link2d ff( 0 ), fn( 0 );
   ff.node = be->n0;
   ff.next = &fn;
   fn.node = be->n1;

   it = fist_front.find( &ff );
   if( it == fist_front.end() )
      lf = *fist_front.begin();
   else
      lf = *it;

   vip1 = lf;
   vim1 = lf->next;

   it = fist_front.find( vim1 );
   if( it == fist_front.end() )
      THROW__X( "mesh2d::create_new_cell: vim1 link not found." );
   fist_front.erase( it );

   it = fist_front.find( vip1 );
   if( it == fist_front.end() )
      THROW__X( "mesh2d::create_new_cell: vip1 link not found." );
   fist_front.erase( it );

   mesh_nodes.push_back( n );
   node2d *nm1 = vim1->node;
   node2d *np1 = vip1->node;

   cell2d *c = cell_alloc( ++cur_cell_id );

   c->face[0].node = np1;
   c->face[1].node = nm1;
   c->face[2].node = n;

   c->face[2].adj = vip1->adj;
   if( vip1->adj )
      vip1->adj->adj = c->face+2;

   l0 = link_alloc( ++cur_link_id );
   l0->node = n;
   l0->next = vim1;
   l0->prev = vip1;
   l0->adj = c->face+0;

   vip1->next = l0;
   vim1->prev = l0;

   vip1->adj = c->face+1;

   circun_circle( c );
   bad_cells.insert( c );
   attach_bad_edges( c );

   fist_front.insert( vip1 );
   fist_front.insert( l0 );
   fist_front.insert( vim1 );
}

bool
mesh2d::is_point_allowed( cell2d *c, node2d *n )
{
   //
   // it is assumed that c contains n
   //
   cell2d_set::iterator itt;
   cell2d_set tested;
   list<face2d*> untested;

   face2d *f, *fa;
   int i;

   if( !c->bad )
      return false;

   if( is_close_to_existing_node( c, n ) )
      return false;
   tested.insert( c );

   for( i = 0 ; i < 3 ; i++ )
   {
      fa = c->face[i].adj;
      if( fa != 0 )
         untested.push_back( fa );
   }

   while( !untested.empty() )
   {
      f = untested.front();
      untested.pop_front();
      c = f->cell;

      itt = tested.find( c );
      if( itt != tested.end() )
         continue;
      else
      {
         if( is_close_to_existing_node( c, n ) )
            return false;
         tested.insert( c );

         if( inside_circuncircle( c, n ) )
         {
            if( !c->bad )
               return false;
            else
            {
               fa = succ_face( f )->adj;
               if( fa )
                  untested.push_back( fa );
               fa = pred_face( f )->adj;
               if( fa )
                  untested.push_back( fa );
            }
         }
      }
   }
   return true;
}

bool
mesh2d::is_close_to_existing_node( cell2d *c, node2d *n )
{
   //
   // it is assumed that c contains n
   //
   double s2, r2;
   node2d *w;
   int i;

   for( i = 0 ; i < 3 ; i++ )
   {
      w  = c->face[i].node;
      s2 = norm_sqr( n->p.x - w->p.x, n->p.y - w->p.y );
      r2 = sqr( dist_factor * 0.5 * ( n->h + w->h ) );
      if( s2 <= r2 )
         return true;
   }
   return false;
}

void
mesh2d::attach_bad_edges( cell2d *c )
{
   edge2d_set_addr::iterator it;
   edge2d e, *a;
   face2d *f;
   int i;

   for( i = 0 ; i < 3 ; i++ )
   {
      f = c->face+i;
      e.n0 = succ_node( f );
      e.n1 = pred_node( f );

      it = frontal_edges.find( &e );
      if( it != frontal_edges.end() )
      {
         a = *it;
         a->adj = f;
      }
   }
}

void
mesh2d::create_frontal_edges()
{
   //
   // A perfectly symmetric algorithm will create a new front, based on the
   // neighbours, without changing the current front. The process will be
   // stopped when all the points do not have neighbours.
   // The current algorithm is enough for our purposes and is more fast.
   //
   edge2d_set_sort::iterator itc, ite;
   edge2d_set_sort frontal_list_1, frontal_list_2;
   edge2d_set_sort *frontal_temp, *frontal_wait;
   edge2d_set_sort neighbours;

   cell2d_set::iterator it;

   cell2d *cp;
   node2d *ne, *nc;
   edge2d *be, *bc;
   face2d *f;

   double s2, r2, ex, ey, eh, ec, eo;
   bool changed;
   int i;

   frontal_edges.clear();
   frontal_temp = &frontal_list_1;
   frontal_wait = &frontal_list_2;

   for( it = bad_cells.begin() ; it != bad_cells.end() ; it++ )
   {
      cp = *it;

      for( i = 0 ; i < 3 ; i++ )
      {
         f = cp->face + i;
         if( is_frontal_face( f ) )
         {
            be = edge_alloc();
            be->new_node = node_alloc( ++cur_node_id );
            be->adj = f;
            be->n0 = succ_node( f );
            be->n1 = pred_node( f );
            be->h  = norm( be->n1->p.x - be->n0->p.x, be->n1->p.y - be->n0->p.y );
            ex = 0.5 * ( be->n1->p.x + be->n0->p.x ) - x_centroid;
            ey = 0.5 * ( be->n1->p.y + be->n0->p.y ) - y_centroid;
            be->theta = atan2pi( ey, ex );
            be->ro = norm( ex, ey );
            new_point( be );

            if( be->ins_cell != 0 && is_point_allowed( be->ins_cell, be->new_node ) )
               frontal_wait->insert( be );
            else
            {
               delete be->new_node;
               delete be;
            }
         }
      }
   }

   do {
      changed = false;
      swap( frontal_wait, frontal_temp );

      while( !frontal_temp->empty() )
      {
         ite = frontal_temp->begin();
         be = *ite;
         ne = be->new_node;
         frontal_temp->erase( ite );

         if( ne == 0 )
         {
            delete be;
            continue;
         }
         //
         // let us collect all the neighbours by
         // checking only the remainder of the set
         //
         neighbours.clear();
         for( itc = frontal_temp->begin() ; itc != frontal_temp->end() ; itc++ )
         {
            bc = *itc;
            nc = bc->new_node;

            if( nc == 0 )
               continue;

            s2 = norm_sqr( ne->p.x - nc->p.x, ne->p.y - nc->p.y );
            r2 = sqr( dist_factor * 0.5 * ( ne->h + nc->h ) );
            if( s2 <= r2 )
               neighbours.insert( bc );
         }
         //
         // evaluate average values
         //
         ex = ne->p.x / ne->h;
         ey = ne->p.y / ne->h;
         eo = 1.0 / ne->h;
         eh = ne->h;
         ec = 1.0;

         if( !neighbours.empty() ) changed = true;

         for( itc = neighbours.begin() ; itc != neighbours.end() ; itc++ )
         {
            bc = *itc;
            nc = bc->new_node;

            if( nc == 0 )
               continue;

            ex += nc->p.x / nc->h;
            ey += nc->p.y / nc->h;
            eo += 1.0 / nc->h;
            eh += nc->h;
            ec += 1.0;
         }
         ne->p.x = ex / eo;
         ne->p.y = ey / eo;
         ne->h = eh / ec;

         be->ins_cell = find_bad_cell( be->ins_cell, ne );
         if( be->ins_cell != 0 &&
             is_point_allowed( be->ins_cell, ne ) )//&&
             //!is_close_to_frontal_face( be->ins_cell, ne ) )
         {
            frontal_wait->insert( be );
            //
            // the neighbours are useless, mark it
            //
            for( itc = neighbours.begin() ; itc != neighbours.end() ; itc++ )
            {
               bc = *itc;
               delete bc->new_node;
               bc->new_node = 0;
            }
         }
         else
         {
            delete ne;
            delete be;
         }
      }

   } while( changed );

   while( !frontal_wait->empty() )
   {
      ite = frontal_wait->begin();
      frontal_edges.insert( *ite );
      frontal_wait->erase( ite );
   }
}

bool
mesh2d::is_implicit( cell2d *c, int type )
{
   double t0x, t0y, t1x, t1y, t2x, t2y;
   double nt0, nt1, nt2, h0, h1, h2;
   double dot01, dot12, dot20, mndot;

   face2d *f;
   int i, flag = 0;

   t0x = c->face[2].node->p.x - c->face[1].node->p.x;
   t0y = c->face[2].node->p.y - c->face[1].node->p.y;
   nt0 = norm( t0x, t0y );

   t1x = c->face[0].node->p.x - c->face[2].node->p.x;
   t1y = c->face[0].node->p.y - c->face[2].node->p.y;
   nt1 = norm( t1x, t1y );

   t2x = c->face[1].node->p.x - c->face[0].node->p.x;
   t2y = c->face[1].node->p.y - c->face[0].node->p.y;
   nt2 = norm( t2x, t2y );

   h0 = impl_factor * 0.5 * ( c->face[1].node->h + c->face[2].node->h );
   h1 = impl_factor * 0.5 * ( c->face[0].node->h + c->face[2].node->h );
   h2 = impl_factor * 0.5 * ( c->face[0].node->h + c->face[1].node->h );

   dot01 = -( t0x * t1x + t0y * t1y );
   dot12 = -( t1x * t2x + t1y * t2y );
   dot20 = -( t2x * t0x + t2y * t0y );
   mndot = min( dot01, min( dot12, dot20 ) );

   for( i = 0 ; i < 3 ; i++ )
   {
      f = c->face[i].adj;
      if( f == 0 || !f->cell->bad )
         flag |= ( 1 << i );
   }

   flag |= ( 1 << ( type + 2 ) );

   switch( flag )
   {
   default:
      return false;
   case 1 +  8:      // type 1
   case 6 + 16:      // type 2
      return( mndot >= 0.0 && nt1 < h1 && nt2 < h2 );
   case 2 +  8:      // type 1
   case 5 + 16:      // type 2
      return( mndot >= 0.0 && nt0 < h0 && nt2 < h2 );
   case 4 +  8:      // type 1
   case 3 + 16:      // type 2
      return( mndot >= 0.0 && nt0 < h0 && nt1 < h1 );
   case 7 + 32:      // type 3
      return true;
   }
}

void
mesh2d::implicit_cells()
{
   cell2d_set::iterator it, itd;
   cell2d *c;
   bool changed = true;
   frontal_edges.clear();
   int type;

   while( changed )
   {
      changed = false;
      for( type = 1 ; type <= 3 ; type++ )
      {
         for( it = bad_cells.begin() ; it != bad_cells.end() ; )
         {
            c = *it;
            itd = it;   // to erase current position if needed
            it++;
            if( is_implicit( c, type ) )
            {
               bad_cells.erase( itd );
               c->bad = false;
               mesh_cells.insert( c );
               changed = true;
            }
         }
      }
   }
}

void
mesh2d::smooth()
{
   cell2d_set::iterator itc;
   cell2d *c;
   node2d_list::iterator itn;
   node2d *nn, *ns, *np;
   int i = 1;

   make_delaunay( &mesh_cells );
   dump_mesh();
   int num = mesh_nodes.size();
   double *nx = new double[num+1];
   double *ny = new double[num+1];
   double *ap = new double[num+1];
   double *sa = new double[num+1];
   //
   // renumber nodes and evaluate node degree
   //
   for( itn = mesh_nodes.begin() ; itn != mesh_nodes.end() ; itn++ )
   {
      nn = *itn;
      nn->degree = 0.0;
      nn->id = i;
      nx[i] = ny[i] = sa[i] = 0.0;
      i++;
   }
   for( itc = mesh_cells.begin() ; itc != mesh_cells.end() ; itc++ )
   {
      c = *itc;
      for( i = 0 ; i < 3 ; i++ )
         c->face[i].node->degree += 1.0;
   }
   //
   //
   //
   for( itn = mesh_nodes.begin() ; itn != mesh_nodes.end() ; itn++ )
   {
      nn = *itn;
      ap[ nn->id ] = max( 6.0, 1.0 + 3.0 * ( nn->degree - 6.0 ) );
   }
   for( itc = mesh_cells.begin() ; itc != mesh_cells.end() ; itc++ )
   {
      c = *itc;
      for( i = 0 ; i < 3 ; i++ )
      {
         nn = c->face[i].node;
         ns = succ_node( &c->face[i] );
         np = pred_node( &c->face[i] );

         if( nn->bc_type != 0 )
            continue;

         sa[ nn->id ] += ap[ ns->id ] + ap[ np->id ];
         nx[ nn->id ] += ap[ ns->id ] * ns->p.x + ap[ np->id ] * np->p.x;
         ny[ nn->id ] += ap[ ns->id ] * ns->p.y + ap[ np->id ] * np->p.y;

      }
   }
   for( itn = mesh_nodes.begin() ; itn != mesh_nodes.end() ; itn++ )
   {
      nn = *itn;

      if( nn->bc_type != 0 )
         continue;

      nn->p.x = nx[ nn->id ] / sa[ nn->id ];
      nn->p.y = ny[ nn->id ] / sa[ nn->id ];
   }

   delete[] nx;
   delete[] ny;
   delete[] ap;
   delete[] sa;

   make_delaunay( &mesh_cells );
}

bool
mesh2d::mesh_generation()
{
   //
   // We assume that we have already a boundary conforming mesh.
   //
   cell2d_set::iterator itb;
   edge2d_set_addr::iterator ita;
   edge2d *be;
   cell2d *ins_cell;
   node2d *new_node;

   int cycle = 1;

   printf( "\nStarting mesh generation\n\n" );

   create_back_mesh();
   dump_back();

   centroide();
   create_frontal_edges();

   while( !frontal_edges.empty() )
   {
      for( ita = frontal_edges.begin() ; ita != frontal_edges.end() ; ita++ )
      {
         be = *ita;

         if( be->adj == 0 ) // this edge has been removed from mesh
         {
            THROW__X( "mesh2d::mesh_generation(): be->adj == 0.\n" );
            delete be->new_node;
            continue;
         }
         //
         // We delete cells each point inserted. Find again...
         //
         new_node = be->new_node;
         ins_cell = find_bad_cell( be->adj->cell, new_node );
         if( ins_cell == 0 )
            THROW__X( "mesh2d::mesh_generation: find_cell failed.\n" );
         //
         // create front for FIST and create a new triangle updating fist_front
         //
         create_fist_front( ins_cell, new_node );
         create_new_cell( be, new_node );
         fist_generation();
      }

      if( make_delaunay( &bad_cells ) ) // should not be needed
        printf( "make_bad_cells_delaunay changed the mesh\n\n" );

      if( !( cycle % 10 ) )
      {
         dump_mesh();
         dump_bad_cells();
      }

      implicit_cells();
      create_frontal_edges();

      printf( "%s   front number = %4i\n", mov_up, cycle++ );
   }
   //
   // move remaining bad_cells to the mesh_cells set
   //
   while( !bad_cells.empty() )
   {
      itb = bad_cells.begin();
      mesh_cells.insert( *itb );
      bad_cells.erase( itb );
   }

   printf( "Ending mesh generation\n\n" );
   return false;
}

void
mesh2d::convert_to_tri6( tri6_xda_interpolator& bi )
{
   int nodes_cur_id, cells_cur_id;

   cell2d* cl = get_first_cell();
   cells_cur_id = 0;
   while( cl )
   {
      cl->id = ++cells_cur_id;
      cl = get_next_cell();
   }

   node2d* nd = get_first_node();
   nodes_cur_id = 0;
   while( nd )
   {
      nd->id = ++nodes_cur_id; // renumbering to load correctly
      nd = get_next_node();
   }

   const int prev_node[3] = { 2, 0, 1 };
   const int next_node[3] = { 1, 2, 0 };
   int i, j, k;

   cl = get_first_cell();
   while( cl )
   {
      for( i = 0; i < 3 ; i++ )
      {
         j = next_node[i];
         k = prev_node[i];

         node2d*  ni = cl->face[i].node;
         node2d*  nj = cl->face[j].node;
         node2d*& nk = cl->face[k].mid_node;

         if( nk == 0 )
         {
            nk = node_alloc( ++nodes_cur_id );
            mid_nodes.push_back( nk );

            if( cl->face[k].adj )
            {
              nk->param = 0.0;
              nk->p.x = 0.5*( ni->p.x + nj->p.x );
              nk->p.y = 0.5*( ni->p.y + nj->p.y );
              cl->face[k].adj->mid_node = nk;
              nk->bc_type = 0;
            }
            else
              bi.interpolate( ni, nj, nk );
         }
      }
      cl = get_next_cell();
   }

}

void
mesh2d::save_tri6_xda( FILE *stream )
{
   int num_cells, num_nodes, num_bc;
   int i;

   cell2d_set&  cs = get_mesh_cells_set();
   node2d_list& nl = get_mesh_nodes_list();
   node2d_list& nm = get_mesh_mid_nodes_list();

   num_cells = cs.size();
   num_nodes = nl.size() + nm.size();
   num_bc = 0;

   for( cell2d_set::iterator itc = cs.begin() ; itc != cs.end() ; itc++ )
   {
      cell2d *cl = *itc;
      for( i = 0; i < 3 ; i++ )
      {
         face2d& f = cl->face[i];
         if( f.adj == 0 )
            num_bc++;
      }
   }

   fprintf( stream, "DEAL 003:003\n" );
   fprintf( stream, "%i\t\t #num elements\n", num_cells );
   fprintf( stream, "%i\t\t #num nodes\n", num_nodes );
   fprintf( stream, "%i\t\t #sum of elem weights\n", 6*num_cells );
   fprintf( stream, "%i\t\t #num of bc\n", num_bc );
   fprintf( stream, "65535\t\t #string size\n" );
   fprintf( stream, "1\t\t #num of elem blocks\n" );
   fprintf( stream, "4\t\t #elem type in the block\n" );
   fprintf( stream, "%i\t\t #num elem in the block\n", num_cells );
   fprintf( stream, "Id string\n" );
   fprintf( stream, "Title string\n" );

   const int vrtx_dump[3] = { 0, 2, 1 };
   const int edge_dump[3] = { 1, 0, 2 };

   for( cell2d_set::iterator itc = cs.begin() ; itc != cs.end() ; itc++ )
   {
      cell2d *cl = *itc;

      for( i = 0; i < 3 ; i++ )
      {
         node2d*  ni = cl->face[ vrtx_dump[i] ].node;
         fprintf( stream, "%i\t", ni->id-1 );
      }
      for( i = 0; i < 3 ; i++ )
      {
         node2d*  ni = cl->face[ edge_dump[i] ].mid_node;
         fprintf( stream, "%i\t", ni->id-1 );
      }
      fprintf( stream, "\n" );
   }

   int max_id = 0;
   for( node2d_list::iterator itl = nl.begin() ; itl != nl.end() ; itl++ )
   {
      node2d *n = *itl;
      fprintf( stream, "%lf\t%lf\t0.0\n", n->p.x, n->p.y );
      max_id = max( max_id, n->id );
   }

   for( node2d_list::iterator itm = nm.begin() ; itm != nm.end() ; itm++ )
   {
      node2d *n = *itm;
      fprintf( stream, "%lf\t%lf\t0.0\n", n->p.x, n->p.y );
      assert( max_id < n->id );
   }

   for( cell2d_set::iterator itc = cs.begin() ; itc != cs.end() ; itc++ )
   {
      cell2d *cl = *itc;

      for( i = 0; i < 3 ; i++ )
      {
         face2d& f = cl->face[ edge_dump[i] ];
         if( f.adj == 0 )
            fprintf( stream, "%i\t%i\t%i\n", cl->id-1, i, f.mid_node->bc_type );
      }
   }

}

//***EOF************************************************************************
