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

#include <vector>
#include <iostream>
#include <fstream>

#include "stopwatch.h"
#include "t_mesh2d.h"
#include "bc2d.h"
#include "getpot.h"

using namespace mesh_2d;

int main( int argc, char* argv[] )
{
   GetPot cmd_ln(argc, argv);

   cout << "argc = " << argc << endl;
   cout << argv[0] << endl;
   cout << argv[1] << endl;
   cout << argv[2] << endl;

   //const bool show = !cmd_ln.search( "-n" );
   const char* ifile = cmd_ln.follow( (char*) 0, 2, "-i","--ifile");
   const char* ofile = cmd_ln.follow( (char*) 0, 2, "-o","--ofile");

   if( ifile == 0 )
   {
     cout << "input file not supplied" << endl;
     cout << "mesh2d_fist -i infile -o outfile" << endl;
     exit(1);
   }

   if( ofile == 0 )
   {
     cout << "output file not supplied" << endl;
     cout << "mesh2d_fist -i infile -o outfile" << endl;
     exit(1);
   }

   mesh2d m;
   stopwatch sw;

   m.set_dumping( false );
   //~ m.set_dump_dir( "./" );

   std::ifstream fs;

   fs.open( ifile );
   if ( !fs )
   {
      std::cerr << "erro na abertura do ficheiro da frente" << std::endl;
      exit(-1);
   }

   double x, y;
   int bc_type = 1, bc_index = 0, bc_surface = 0, n_pts, n_frts;

   fs >> n_frts;

   for( int f = 0; f < n_frts; f++ )
   {
     fs >> n_pts;

     m.begin_front();
     for( int j = 0 ; j < n_pts ; j++ )
     {
        fs >> x >> y >> bc_type;

        if( j != n_pts-1 )
           m.add_to_front( x, y, 0.0, bc_type, bc_index, bc_surface );

        if( !fs )
        {
           std::cout << "erro na leitura do ficheiro da frente (" << f << ")" << std::endl;
           exit(-1);
        }
     }
     m.end_front();
   }
   fs.close();

   sw.start();
   //~ m.dump_links();
   m.fist_generation();
   //~ m.dump_mesh();
   printf( "elapsed time = %.2fs\n", sw.stop() );

   sw.start();
   m.mesh_generation();
   printf( "elapsed time = %.2fs\n", sw.stop() );

   m.testing_mesh();

   //~ m.dump_mesh();
   //~ m.dump_back();
   //~ m.dump_bad_cells();

   m.smooth();
   m.smooth();
   m.smooth();
   m.smooth();
   m.smooth();
   //~ m.dump_mesh();

   FILE* stream = fopen( ofile, "w+" );

   //m.save( stream );
   m.save_gmsh( stream );
   fclose( stream );

   printf( "SUCCESS!\n\n" );
   return 0;
}
//***EOF************************************************************************
