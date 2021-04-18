/*
 * levmarq.c
 *
 * This file contains an implementation of the Levenberg-Marquardt algorithm
 * for solving least-squares problems, together with some supporting routines
 * for Cholesky decomposition and inversion.  No attempt has been made at
 * optimization.  In particular, memory use in the matrix routines could be
 * cut in half with a little effort (and some loss of clarity).
 *
 * It is assumed that the compiler supports variable-length arrays as
 * specified by the C99 standard.
 *
 * Ron Babich, May 2008
 *
 * levmarq.c, levmarq.h, and examples are provided under the MIT license.

   Copyright (c) 2008-2016 Ron Babich

   Modified for C++/std::library: Bernd Wokoeck
   see https://stackoverflow.com/questions/16004668/c-allocating-a-matrix-in-a-function/27366086#27366086
 */

#include <vector>

class LevMarq
{
public:
    struct point_t
    {
        point_t() { x = 0.0; y = 0.0; }
        point_t( double ax, double ay ) { x = ax; y = ay; }
        double x;
        double y;
    };

    typedef std::vector<point_t> measurePoints;

    // cholesky
    typedef std::vector< std::vector<double> > Matrix;
    typedef std::vector<double> Vector;

protected:
    typedef struct {
      int verbose;
      int max_it;
      double init_lambda;
      double up_factor;
      double down_factor;
      double target_derr;
      int final_it;
      double final_err;
      double final_derr;
    } LMstat;

    void levmarq_init(LMstat *lmstat);

    int levmarq( Vector & par, const measurePoints & points, double *dysq,
            double (*func)(Vector&, double, void *),
            void (*grad)(Vector&, Vector&, double, void *),
            void *fdata, LMstat *lmstat);

    double error_func(Vector& par, const measurePoints &points, double *dysq,
          double (*func)(Vector&, double, void *), void *fdata);

    void solve_axb_cholesky(int n, Matrix& l, Vector& x, Vector & b);

    int cholesky_decomp(int n, Matrix& l, Matrix& a );
};
