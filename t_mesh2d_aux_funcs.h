/***************************************************************************
                            mesh generation code
                            --------------------
    First version   [ 001 ] : Oct 17, 2001, 20:02
    Current version [ 088 ] : Jan 05, 2002, 14:15
    copyright               : (C) 2001 by Joao Carlos de Campos Henriques
    email                   : jcch@popsrv.ist.utl.pt
 ***************************************************************************/

inline int
isgn( double a )
{
   return( a > 0.0 ? +1 : ( a < 0.0 ? -1 : 0 ) );
}

inline double
sqr( double x )
{
   return x*x;
}

inline double
atan2pi( double y, double x )
{
   double t = atan2( y, x );
   return( t < 0.0 ? t + 2.0 * M_PI : t );
}

inline double
norm_sqr( double x, double y )
{
   return ( x*x + y*y );
}

inline double
norm( double x, double y )
{
   return sqrt( x*x + y*y );
}

inline void
circun_circle( double X0, double Y0, double X1, double Y1, double X2, double Y2,
               double &xc, double &yc, double &rc )
{
#if 1
   double X01 = X0 - X1;
   double X20 = X2 - X0;

   double Y01 = Y0 - Y1;
   double Y12 = Y1 - Y2;
   double Y20 = Y2 - Y0;
   //
   // evaluated with Mathematica 4.0
   //
   xc = (sqr(X2)*Y01 + sqr(X1)*Y20 + Y12*(sqr(X0) - Y01*Y20))/(2.0*(X2*Y01 + X0*Y12 + X1*Y20));

   yc = -(2.0*X20*(sqr(X0) - sqr(X1) + sqr(Y0) - sqr(Y1)) +
          2.0*X01*(sqr(X0) - sqr(X2) + sqr(Y0) - sqr(Y2)))/(4.0*(-(X2*Y01) - X0*Y12 - X1*Y20));

   rc = norm( X0 - xc, Y0 - yc );
#else
   double a,b,c,d,e,f,g,h,i,j,k;

   a = X0 - X1;
   b = X0 + X1;
   c = Y0 - Y1;
   d = Y0 + Y1;
   e = 0.5 * (a*b + c*d);
   f = X1 - X2;
   g = X1 + X2;
   h = Y1 - Y2;
   i = Y1 + Y2;
   j = 0.5 * (f*g + h*i );
   k = a*h-c*f;

   if ( k == 0.0 || ( a == 0.0 && f == 0.0 ) )
      THROW__X( "circun_circle: colinear points." );

   yc = (a*j-e*f)/k;
   xc = ( fabs(a) > fabs(f) ? (e-c*yc)/a : (j-h*yc)/f );
   rc = norm( X1 - xc, Y1 - yc );
#endif
}

inline int
det( node2d *A, node2d *B, node2d* C )
{
   double Dx  = B->p.x - A->p.x;
   double Dy  = A->p.y - B->p.y;
   double Dxy = A->p.x * B->p.y - A->p.y * B->p.x;
   double dt  = C->p.x * Dy + C->p.y * Dx + Dxy;

   if( fabs( dt ) <= mesh_eps )
      return 0;
   else
      return isgn( dt );
}

inline void
circun_circle( cell2d *cl )
{
   double p0x = cl->face[0].node->p.x;
   double p0y = cl->face[0].node->p.y;

   double p1x = cl->face[1].node->p.x;
   double p1y = cl->face[1].node->p.y;

   double p2x = cl->face[2].node->p.x;
   double p2y = cl->face[2].node->p.y;

   circun_circle( p0x, p0y, p1x, p1y, p2x, p2y, cl->xc, cl->yc, cl->rc );
}

inline node2d*
succ_node( face2d *f )
{
   int id = f->id;
   return f->cell->face[ succ[id] ].node;
}

inline node2d*
pred_node( face2d *f )
{
   int id = f->id;
   return f->cell->face[ pred[id] ].node;
}

inline face2d*
succ_face( face2d *f )
{
   return f->cell->face + succ[ f->id ];
}

inline face2d*
pred_face( face2d *f )
{
   return f->cell->face + pred[ f->id ];
}

inline double
cell_area( cell2d* c )
{
   double x_1 = c->face[0].node->p.x;
   double y_1 = c->face[0].node->p.y;

   double x_2 = c->face[1].node->p.x;
   double y_2 = c->face[1].node->p.y;

   double x_3 = c->face[2].node->p.x;
   double y_3 = c->face[2].node->p.y;

   return -0.5*( x_2*y_3 + x_1*y_2 + x_3*y_1 - y_1*x_2 - y_2*x_3 - y_3*x_1 );
}

inline bool
segs_intersection( node2d *p0, node2d *p1, node2d *q0, node2d *q1 )
{
   double d1x = p1->p.x - p0->p.x;
   double d1y = p1->p.y - p0->p.y;

   double d2x = q1->p.x - q0->p.x;
   double d2y = q1->p.y - q0->p.y;

   double d0x = q0->p.x - p0->p.x;
   double d0y = q0->p.y - p0->p.y;

   double den = -d1x*d2y + d2x*d1y;

   if( den == 0.0 )
      return false;

   double t1 = ( -d2y*d0x + d2x*d0y ) / den;
   double t2 = ( -d1y*d0x + d1x*d0y ) / den;

   return ( 0.0 <= t1 && t1 <= 1.0 && 0.0 <= t2 && t2 <= 1.0 );
}

inline void
normal_to_face( node2d* V[2], p2d* n )
{
   n->x = + V[1]->p.y - V[0]->p.y;
   n->y = - V[1]->p.x + V[0]->p.x;
}

inline void
cell_normals( cell2d* tc )
{
   node2d* nd[2];
   int f, t;
   for( f = 0 ; f < 3 ; f++ )
   {
      for( t = 0 ; t < 2 ; t++ )
         nd[t] = tc->face[ FAC2D[f][t] ].node;
      normal_to_face( nd, &tc->face[f].normal );
   }
}

inline int get_cell_type( cell2d *tc )
{
  union ctype {
    struct {    
      int inv_face_0 : 1;
      int inv_face_1 : 1;
      int inv_face_2 : 1;
    } bits;
    int type;
  };

  ctype tt;  
  
  assert( tc->face[0].orientation != 0 );
  assert( tc->face[1].orientation != 0 );
  assert( tc->face[2].orientation != 0 );

  tt.type = 0;
  tt.bits.inv_face_0 = ( tc->face[0].orientation == -1 );
  tt.bits.inv_face_1 = ( tc->face[1].orientation == -1 );
  tt.bits.inv_face_2 = ( tc->face[2].orientation == -1 );
  
  assert( tt.type > 0 &&  tt.type < 7 );

  return tt.type;
}
//***EOF************************************************************************
