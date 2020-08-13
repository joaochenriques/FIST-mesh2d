

/***************************************************************************
                          common.h  -  description
                             -------------------
    begin                : Tue Feb 19, 2001
    copyright            : (C) 2001 by Joao Carlos de Campos Henriques
    email                : jcch@hidro1.ist.utl.pt
 ***************************************************************************/
#ifndef COMMON_H
#define COMMON_H

#include <string.h>
#include <stdlib.h>
#include <math.h>

#define mem_zero(a,b)     memset( a, 0, sizeof(b) )

struct __X
{
   const char* error_msg;
   const char* at_file;
   const int   at_line;

   __X( const char *s, const char *f, const int l ) : error_msg(s), at_file(f), at_line(l) {}
};

#define THROW__X( a )            \
   {                             \
      static const char __s[] = a;     \
      static __X __x(__s,__FILE__,__LINE__);       \
      printf( "\n\nABORTING PROGRAM at file: %s in line %i\n\n%s\n\n", \
              __FILE__, __LINE__, a ); \
      exit(-1);                  \
   }

#define ASSERT__X( a )                                \
if( !(a) )                                            \
   {                                                     \
      static const char __s[] = "assertion ( " #a " ) failed.";\
      static __X __x(__s,__FILE__,__LINE__);          \
      printf( "\n\nABORTING PROGRAM at file: %s in line %i\n\n%s\n\n", \
              __FILE__, __LINE__, __s ); \
      exit(-1);                  \
   }

template<class T>
inline void sink(T&) {};

template<class type>
inline type tsqr( type a )
{
   return a * a;
}

template<class type>
inline type tmax( type a, type b )
{
   return( a > b ? a : b );
}

template<class type>
inline type tmin( type a, type b )
{
   return( a < b ? a : b );
}

template<class T>
class NRVec {
protected:
  int nn;
  T *v;
public:
  NRVec() : nn(0), v(0) {}

  explicit NRVec( int n ) : nn(n), v(new T[n]) {}

  NRVec( const T& a, int n ) : nn(n), v(new T[n])
  {
    int i;
    for( i = 0 ; i < nn ; i++ )
      v[i] = a;
  }

  NRVec( const T* a, int n ) : nn(n), v(new T[n])
  {
    int i;
    for( i = 0 ; i < nn ; i++ )
      v[i] = *a++;
  }

  NRVec( const NRVec &rhs ) : nn(rhs.nn), v(new T[rhs.nn])
  {
    int i;
    for( i = 0 ; i < nn ; i++ )
      v[i] = rhs[i];
  }

  void resize( int n )
  {
    if( v != 0 )
      delete[] v;
    nn = n;
    v  = new T[nn];
  }

  NRVec& operator = ( const NRVec &rhs )
  {
    if( this != &rhs )
    {
      if( nn != rhs.nn )
        resize( rhs.nn );
      int i;
      for( i = 0 ; i < nn ; i++ )
        v[i] = rhs[i];
    }
    return *this;
  }

  NRVec& operator = ( const T& a )
  {
    int i;
    for( i = 0 ; i < nn ; i++ )
      v[i] = a;
    return *this;
  }

  inline T& operator[] ( int i )
  {
    return v[i];
  }

  inline const T& operator[] ( const int i ) const
  {
    return v[i];
  }

  inline int size() const
  {
    return nn;
  }

  ~NRVec()
  {
    if( v != 0 )
      delete[] v;
  }
};

template<class T>
class NRMat {
protected:
  int nn;
  int mm;
  T **v;
public:
  NRMat() : nn(0), mm(0), v(0) {}

  NRMat( int n, int m ) : nn(n), mm(m), v(new T*[n])
  {
    v[0] = new T[m*n];
    for( int i = 1 ; i < n ; i++ )
      v[i] = v[i-1] + m;
  }

  NRMat( const T& a, int n, int m ) : nn(n), mm(m), v(new T*[n] )
  {
    int i, j;
    v[0] = new T[m*n];
    for( i = 1 ; i < n ; i++ )
      v[i] = v[i-1] + m;
    for( i = 0 ; i < n ; i++ )
      for( j = 0 ; j < m ; j++ )
        v[i][j] = a;
  }

  NRMat( const T* a, int n, int m ) : nn(n), mm(m), v(new T*[n] )
  {
    int i, j;
    v[0] = new T[m*n];
    for( i = 1 ; i < n ; i++ )
      v[i] = v[i-1] + m;
    for( i = 0 ; i < n ; i++ )
      for( j = 0 ; j < m ; j++ )
        v[i][j] = *a++;
  }

  NRMat( const NRMat &rhs ) : nn(rhs.nn), mm(rhs.mm), v(new T*[rhs.nn])
  {
    int i, j;
    v[0] = new T[nn*mm];
    for( i = 1 ; i < nn ; i++ )
      v[i] = v[i-1] + mm;

    for( i = 0 ; i < nn ; i++ )
      for( j = 0 ; j < mm ; j++ )
        v[i][j] = rhs[i][j];
  }

  void resize( int n, int m )
  {
    if( v != 0 )
    {
      delete[] v[0];
      delete[] v;
    }
    nn = n;
    mm = m;
    v = new T*[nn];
    v[0] = new T[mm*nn];
    for( int i = 1 ; i < nn ; i++ )
      v[i] = v[i-1] + mm;
  }

  NRMat& operator = ( const NRMat &rhs )
  {
    if( this != &rhs )
    {
      int i, j;
      if( nn != rhs.nn || mm != rhs.mm )
        resize( rhs.nn, rhs.mm );
      for( i = 0 ; i < nn ; i++ )
        for( j = 0 ; j < mm ; j++ )
          v[i][j] = rhs[i][j];
    }
    return *this;
  }

  inline NRMat& operator = ( const T& a )
  {
    int i, j;
    for( i = 0 ; i < nn ; i++ )
      for( j = 0 ; j < mm ; j++ )
        v[i][j] = a;
    return *this;
  }

  inline T* operator [] ( const int i )
  {
    return v[i];
  }

  inline const T* operator [] ( const int i )  const
  {
    return v[i];
  }

  inline int nrows() const
  {
    return nn;
  }

  inline int ncols() const
  {
    return mm;
  }

  ~NRMat()
  {
     if( v != 0 )
     {
       delete[] v[0];
       delete[] v;
     }
  }
};

typedef NRVec<int> Vec_I;
typedef NRVec<double> Vec_DP;
typedef NRMat<double> Mat_DP;
typedef NRVec<Vec_DP> VVec_DP;

/*class ilu_solve
{
 protected:
  double tiny;
  Vec_I  indx;
  int    N;
 public:
  ilu_solve( int n ) : tiny(1E-20), indx(n), N(n) {}

  ~ilu_solve() {}

  bool decomposition( Mat_DP &a )
    {
      int i, imax, j, k;
      double big, dum, sum, temp;
      Vec_DP vv( N );

      assert( a.nrows() == a.ncols() );
      assert( N == a.ncols() );

      imax = -1;

      for( i = 0 ; i < N ; i++ )
      {
        big = 0.0;
        for( j = 0 ; j < N ; j++ )
        {
          if( ( temp = fabs(a[i][j]) ) > big )
          {
            big = temp;
          }
        }
        if( big == 0.0 )
        {
          return true;
        }
        vv[i] = 1.0 / big;
      }
      for( j = 0 ; j < N ; j++ )
      {
        for( i = 0 ; i < j ; i++ )
        {
          sum = a[i][j];
          for( k = 0 ; k < i ; k++ )
          {
            sum -= a[i][k] * a[k][j];
          }
          a[i][j] = sum;
        }
        big = 0.0;
        for( i = j ; i < N ; i++ )
        {
          sum = a[i][j];
          for( k = 0 ; k < j ; k++ )
          {
            sum -= a[i][k] * a[k][j];
          }
          a[i][j] = sum;
          if( ( dum = vv[i] * fabs( sum ) ) >= big )
          {
            big = dum;
            imax = i;
          }
        }
        if( j != imax )
        {
          for( k = 0 ; k < N ; k++ )
          {
            dum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = dum;
          }
          vv[imax] = vv[j];
        }
        indx[j] = imax;
        if( a[j][j] == 0.0 )
        {
          a[j][j] = tiny;
        }
        if( j != N - 1 )
        {
          dum = 1.0 / ( a[j][j] );
          for(i = j + 1 ; i < N ; i++ )
          {
            a[i][j] *= dum;
          }
        }
      }
      return false;
    }

  void back_substitution( Mat_DP &a, Vec_DP &b )
    {
      int i, ii = -1, ip, j;
      double sum;

      assert( a.nrows() == a.ncols() );
      assert( N == a.ncols() );
      assert( N == b.size() );

      for( i = 0; i < N; i++ )
      {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if( ii > -1 )
        {
          for( j = ii ; j <= i - 1 ; j++ )
          {
            sum -= a[i][j] * b[j];
          }
        } else if( sum ) {
          ii = i;
        }
        b[i] = sum;
      }
      for( i = N - 1 ; i >= 0 ; i-- )
      {
        sum = b[i];
        for( j = i + 1 ; j < N ; j++ )
        {
          sum -= a[i][j] * b[j];
        }
        b[i] = sum / a[i][i];
      }
    }
};*/

#endif

/*
* end-of-file
*/
