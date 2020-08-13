#ifndef STOPWATCH_H
#define STOPWATCH_H

// for clock() and CLOCKS_PER_SEC
#include <time.h>

/*  Simple stopwatch object:

        void    start()     : start timing
        double  stop()      : stop timing
        void    reset()     : set elapsed time to 0.0
        double  read()      : read elapsed time (in seconds)
*/

class stopwatch
{
private:
   bool   running;
   time_t last_time;
   double total;

public:
   stopwatch() : running(0), last_time(0), total(0.0) {}

   void reset()
   {
      running = 0;
      last_time = 0;
      total=0.0;
   }

   void start()
   {
      if( !running )
      {
         last_time = time(0);
         running = true;
      }
   }

   double stop()
   {
      if( running )
      {
         total += difftime( time(0), last_time );
         running = false;
      }
      return total;
   }

   double read()
   {
      if( running )
      {
         total += difftime( time(0), last_time );
         last_time = time(0);
      }
      return total;
   }

};

#endif
