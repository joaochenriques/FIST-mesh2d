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

#define PATH_MAX 512

#include "t_mesh2d.h"

using namespace mesh_2d;

bool
mesh2d::check_dump_dir()
{
  if( dump_dir == 0 )
  {
    printf("mesh2d::dump_dir not defined.\nUnable to dump\n" );
    return false;
  } else
    return true;
}

void
mesh2d::dump_mesh()
{
   if( !dumping_enable ) return;

   if( !check_dump_dir() ) return;

   char filename[PATH_MAX];
   strcpy( filename, dump_dir );
   strcat( filename, "dump_mesh.txt" );

   FILE* res = fopen( filename, "w" );
   check_stream( res, filename );

   cell2d_set::iterator it;
   cell2d *c;

   for( it = mesh_cells.begin() ; it != mesh_cells.end() ; it++ )
   {
      c = *it;
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[1].node->p.x, c->face[1].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[2].node->p.x, c->face[2].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "\n" );
   }

   for( it = bad_cells.begin() ; it != bad_cells.end() ; it++ )
   {
      c = *it;
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[1].node->p.x, c->face[1].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[2].node->p.x, c->face[2].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
#if 0
      double xm = ( c->face[0].node->x + c->face[1].node->p.x + c->face[2].node->p.x ) / 3.0;
      double ym = ( c->face[0].node->y + c->face[1].node->p.y + c->face[2].node->p.y ) / 3.0;
      fprintf( res, "%f\t%f\n", xm, ym );
#endif
      fprintf( res, "\n" );
   }

   fflush( res );
   fclose( res );

   printf( "dump_mesh\n" );
   printf( "number of cells = %i\n", mesh_cells.size() + bad_cells.size() );
   printf( "number of nodes = %i\n", mesh_nodes.size() );
   printf( "number of convex links = %i\n", convex.size() );
   printf( "number of reflex links = %i\n", reflex.size() );
   printf( "number of waiting links = %i\n", waiting.size() );
   printf( "number of links = %i\n", fist_front.size() );
}

void
mesh2d::dump_bad_cells()
{
   if( !dumping_enable ) return;

   if( !check_dump_dir() ) return;

   char filename[PATH_MAX];
   strcpy( filename, dump_dir );
   strcat( filename, "dump_bad.txt" );

   FILE* res = fopen( filename, "w" );
   check_stream( res, filename );

   cell2d_set::iterator it;
   cell2d *c;

   for( it = bad_cells.begin() ; it != bad_cells.end() ; it++ )
   {
      c = *it;
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[1].node->p.x, c->face[1].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[2].node->p.x, c->face[2].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "\n" );
   }

   fflush( res );
   fclose( res );
}

void
mesh2d::dump_back()
{
   if( !dumping_enable ) return;

   if( !check_dump_dir() ) return;

   char filename[PATH_MAX];
   strcpy( filename, dump_dir );
   strcat( filename, "dump_back.txt" );

   FILE* res = fopen( filename, "w" );
   check_stream( res, filename );

   cell2d_set::iterator it;
   cell2d *c;

   for( it = back_mesh.begin(); it != back_mesh.end() ; it++ )
   {
      c = *it;
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[1].node->p.x, c->face[1].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[2].node->p.x, c->face[2].node->p.y );
      fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
      fprintf( res, "\n" );
   }
   fflush( res );
   fclose( res );
}

void
mesh2d::dump_links()
{
   if( !dumping_enable ) return;

   if( !check_dump_dir() ) return;

   char filename[PATH_MAX];
   strcpy( filename, dump_dir );
   strcat( filename, "dump_links.txt" );

   FILE* res = fopen( filename, "w" );
   check_stream( res, filename );

   link2d_set_addr::iterator it;
   link2d *l;

   for(it = fist_front.begin() ; it != fist_front.end() ; it++ )
   {
      l = *it;
      fprintf( res, "%f\t%f\n", l->node->p.x, l->node->p.y );
      fprintf( res, "%f\t%f\n", l->next->node->p.x, l->next->node->p.y );
      fprintf( res, "\n" );
   }
   fflush( res );
   fclose( res );
}

void
mesh2d::dump_frontal_edges()
{
   if( !dumping_enable ) return;

   if( !check_dump_dir() ) return;

   char filename[PATH_MAX];
   strcpy( filename, dump_dir );
   strcat( filename, "dump_front.txt" );

   FILE* res = fopen( filename, "w" );
   check_stream( res, filename );

   edge2d_set_addr::iterator it;
   edge2d *l;

   for( it = frontal_edges.begin() ; it != frontal_edges.end() ; it++ )
   {
      l = *it;
      if( l->adj != 0 )
      {
         fprintf( res, "%f\t%f\n", l->n0->p.x, l->n0->p.y );
         fprintf( res, "%f\t%f\n", l->n1->p.x, l->n1->p.y );
         fprintf( res, "%f\t%f\n", l->new_node->p.x, l->new_node->p.y );
         fprintf( res, "\n" );
      }
   }
   fflush( res );
   fclose( res );
}

void
mesh2d::dump_cell( cell2d *c, cell2d *s, node2d *n )
{
   if( !dumping_enable ) return;

   if( !check_dump_dir() ) return;

   char filename[PATH_MAX];
   strcpy( filename, dump_dir );
   strcat( filename, "dump_cell.txt" );

   FILE* res = fopen( filename, "w" );
   check_stream( res, filename );

   cell2d *a;
   double xm = 0.0, ym = 0.0;
   int i, j;

   for( i = 0 ; i < 3 ; i++ )
   {
      fprintf( res, "%f\t%f\n", c->face[i].node->p.x, c->face[i].node->p.y );
      xm += c->face[i].node->p.x / 3.0;
      ym += c->face[i].node->p.y / 3.0;
   }
   fprintf( res, "%f\t%f\n", c->face[0].node->p.x, c->face[0].node->p.y );
   fprintf( res, "%f\t%f\n", xm, ym );
   fprintf( res, "\n" );

   for( j = 0 ; j < 3 ; j++ )
      if( c->face[j].adj )
      {
         a = c->face[j].adj->cell;
         for( i = 0 ; i < 3 ; i++ )
            fprintf( res, "%f\t%f\n", a->face[i].node->p.x, a->face[i].node->p.y );
         fprintf( res, "%f\t%f\n", a->face[0].node->p.x, a->face[0].node->p.y );
         fprintf( res, "\n" );
      }

   fprintf( res, "%f\t%f\n", c->xc, c->yc );
   fprintf( res, "%f\t%f\n", n->p.x, n->p.y );
   fprintf( res, "%f\t%f\n", xm, ym );
   fprintf( res, "\n" );

   for( i = 0 ; i < 3 ; i++ )
      fprintf( res, "%f\t%f\n", s->face[i].node->p.x, s->face[i].node->p.y );
   fprintf( res, "%f\t%f\n", s->face[0].node->p.x, s->face[0].node->p.y );

   fflush( res );
   fclose( res );
}

void
mesh2d::testing_mesh()
{
   cell2d_set::iterator it = mesh_cells.begin();
   cell2d *cl;
   int i;

   printf( "\nmesh2d::testing_mesh: testing mesh\n" );
   while( it != mesh_cells.end() )
   {
      cl = *it;

      assert( cell_area( cl ) > 0.0 );

      for( i = 0 ; i < 3 ; i++ )
      {
         if( cl->face[i].adj != 0 && cl->face[i].adj->adj != (cl->face+i) )
            printf( "mesh2d::testing_mesh: invalid link\n" );
      }
      it++;
   }
   printf( "mesh2d::testing_mesh: done\n\n" );
}

bool
mesh2d::is_delaunay()
{
   cell2d_set::iterator it = mesh_cells.begin();
   cell2d *c;
   face2d *f;
   node2d *n;
   double  s;
   int i, cnt = 0;

   while( it != mesh_cells.end() )
   {
      c = *it;
      for( i = 0 ; i < 3 ; i++ )
      {
         f = c->face[i].adj;
         if( f != 0 )
         {
            n = f->node;
            s  = norm( n->p.x - c->xc, n->p.y - c->yc );
            if( s <= c->rc * ( 1.0 + 1E-8 ) )
               cnt++;
         }
      }
      it++;
   }
   printf( "Non delaunay mesh. %i bad mesh cells.\n", cnt );
   return( cnt == 0 );
}

bool
mesh2d::green_sibson( face2d *f1 )
{
   node2d *tn1[3];
   node2d *tn2[3];

   face2d *tf1[3];
   face2d *tf2[3];

   cell2d *c1, *c2;
   face2d *f2;
   node2d *n2;
   double s, a, b;
   int i;

   if( f1->adj != 0 )
   {
      c1 = f1->cell;
      f2 = f1->adj;
      n2 = f2->node;

      s  = norm( c1->xc - n2->p.x, c1->yc - n2->p.y );

      if( s <= c1->rc * ( 1.0 - 1E-8 ) )
      {
         // non-Delaunay
         c2 = f2->cell;

         tn1[0] = succ_node( f1 );
         tn1[1] = n2;
         tn1[2] = f1->node;

         tn2[0] = pred_node( f1 );
         tn2[1] = f1->node;
         tn2[2] = n2;

         a = det( tn1[2], tn1[0], tn1[1] );
         b = det( tn2[1], tn2[0], tn2[2] );

         if( a >= 0 || b <= 0 )
            return false;

         tf1[0] = c2->face+0;
         tf1[1] = pred_face( f1 )->adj;
         tf1[2] = succ_face( f2 )->adj;

         tf2[0] = c1->face+0;
         tf2[1] = pred_face( f2 )->adj;
         tf2[2] = succ_face( f1 )->adj;

         for( i = 0 ; i < 3 ; i++ )
         {
            c1->face[i].node = tn1[i];
            c2->face[i].node = tn2[i];

            attach_faces( c1->face+i, tf1[i] );
            attach_faces( c2->face+i, tf2[i] );
         }
         circun_circle( c1 );
         circun_circle( c2 );
         if( tf1[2] )
            green_sibson( tf1[2] );
         if( tf2[1] )
            green_sibson( tf2[1] );
         return true;
      }
      else
         return false;
   }
   return false;
}

bool
mesh2d::make_delaunay( cell2d_set* s )
{
   cell2d_set::iterator it;
   cell2d *c;
   bool changed = true;
   int i, cnt = 0;

   if( s->empty() )
      THROW__X( "mesh2d::make_delaunay called with empty set.\n" );

   while( changed )
   {
      changed = false;

      for( it = s->begin() ; it != s->end() ; it++ )
      {
         c = *it;
         for( i = 0 ; i < 3 ; i++ )
            changed |= green_sibson( c->face+i );
      }
      cnt++;
   }
   return( cnt == 1 );
}
//***EOF************************************************************************
