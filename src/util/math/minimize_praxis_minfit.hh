// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//

// $Id: minimize_praxis_minfit.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "math_util_exception.h"
#include "minimize_praxis_minfit.h"

#include <cmath>


namespace PRAXIS {

// C...AN IMPROVED VERSION OF MINFIT (SEE GOLUB AND REINSCH, 1969)
// C   RESTRICTED TO M=N,P=0.
// C   THE SINGULAR VALUES OF THE ARRAY AB ARE RETURNED IN Q AND AB IS
// C   OVERWRITTEN WITH THE ORTHOGONAL MATRIX V SUCH THAT U.DIAG(Q) = AB.V,
// C   WHERE U IS ANOTHER ORTHOGONAL MATRIX.
//
template <typename FloatType>
void
minfit(const int n,
       FloatType eps,
       const FloatType tol,
       FloatType* ab,
       FloatType* q,
       FloatType* minfit_ws){

  int l;
  FloatType c, f, h, s, y, z;
  FloatType* e(minfit_ws);

  /* householder's reduction to bidiagonal form */
  FloatType x = 0.;
  FloatType g = 0.;
  for (int i=0; i<n; i++) {
    e[i] = g;
    s = 0.0;
    l = i+1;
    for (int j=i; j<n; j++)
      s += ab[j+i*n] * ab[j+i*n];
    if (s < tol) {
      g = 0.0;
    } else {
      f = ab[i+i*n];
      if (f < 0.0)
        g = std::sqrt(s);
      else
        g = -std::sqrt(s);
      h = f*g - s;
      ab[i+i*n] = f - g;
      for (int j=l; j<n; j++) {
	      f = 0.0;
	      for (int k=i; k<n; k++)
          f += ab[k+i*n] * ab[k+j*n];
	      f /= h;
	      for (int k=i; k<n; k++)
          ab[k+j*n] += f * ab[k+i*n];
      }
    }
    q[i] = g;
    s = 0.0;
    if (i < n)
      for (int j=l; j<n; j++)
	      s += ab[i+j*n] * ab[i+j*n];
    if (s < tol) {
      g = 0.0;
    } else {
      f = ab[i+(i+1)*n];
      if (f < 0.0)
        g = std::sqrt(s);
      else
        g = - std::sqrt(s);
      h = f*g - s;
      ab[i+(i+1)*n] = f - g;
      for(int j=l; j<n; j++)
	      e[j] = ab[i+j*n]/h;
      for (int j=l; j<n; j++) {
	      s = 0;
	      for(int k=l; k<n; k++) s += ab[j+k*n]*ab[i+k*n];
	      for(int k=l; k<n; k++) ab[j+k*n] += s * e[k];
      }
    }
    y = std::fabs(q[i]) + std::fabs(e[i]);
    if (y > x) x = y;
  }
  /* accumulation of right hand transformations */
  for (int i=n-1; i >= 0; i--) {
    if (g != 0.0) {
      h = ab[i+(i+1)*n]*g;
      for (int j=l; j<n; j++) ab[j+i*n] = ab[i+j*n] / h;
      for (int j=l; j<n; j++) {
        s = 0.0;
	      for (int k=l; k<n; k++) s += ab[i+k*n] * ab[k+j*n];
	      for (int k=l; k<n; k++) ab[k+j*n] += s * ab[k+i*n];
      }
    }
    for (int j=l; j<n; j++)
      ab[i+j*n] = ab[j+i*n] = 0.0;
    ab[i+i*n] = 1.0;
    g = e[i];
    l = i;
  }
  /* diagonalization to bidiagonal form */
  eps *= x;
  for (int k=n-1; k>= 0; k--) {
    int kt = 0;
TestFsplitting:
    if (++kt > 30) {
      e[k] = 0.0;
      throw math_util_exception("+++ qr failed");
    }
    for (int l2=k; l2>=0; l2--) {
      l = l2;
      if (std::fabs(e[l]) <= eps)
	      goto TestFconvergence;
      if (std::fabs(q[l-1]) <= eps)
        break;	/* goto Cancellation; */
    }
    //C...CANCELLATION OF E(L) IF L>1...
//Cancellation:
    c = 0.0;
    s = 1.0;
    for (int i=l; i<=k; i++) {
      f = s * e[i];
      e[i] *= c;
      if (std::fabs(f) <= eps)
	      goto TestFconvergence;
      g = q[i];
      if (std::fabs(f) < std::fabs(g)) {
	      FloatType fg = f/g;
	      h = std::fabs(g)*std::sqrt(1.0+fg*fg);
      } else {
	      FloatType gf = g/f;
	      h = (f!=0.0 ? fabs(f)*std::sqrt(1.0+gf*gf) : 0.0);
      }
      q[i] = h;
      if (h == 0.0) { h = 1.0; g = 1.0; }
      c = g/h; s = -f/h;
    }
TestFconvergence:
    z = q[k];
    if (l == k)
      goto Convergence;
    /* shift from bottom 2x2 minor */
    x = q[l];
    y = q[k-l];
    g = e[k-1];
    h = e[k];
    f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h*y);
    g = std::sqrt(f*f+1.0);
    if (f <= 0.0)
      f = ((x-z)*(x+z) + h*(y/(f-g)-h))/x;
    else
      f = ((x-z)*(x+z) + h*(y/(f+g)-h))/x;
    /* next qr transformation */
    s = c = 1.0;
    for (int i=l+1; i<=k; i++) {
      g = e[i];
      y = q[i];
      h = s*g;
      g *= c;
      if (std::fabs(f) < std::fabs(h)) {
	      FloatType fh = f/h;
	      z = fabs(h) * sqrt(1.0 + fh*fh);
      } else {
	      FloatType hf = h/f;
	      z = (f!=0.0 ? std::fabs(f)*std::sqrt(1.0+hf*hf) : 0.0);
      }
      e[i-1] = z;
      if (z == 0.0)
 	      f = z = 1.0;
      c = f/z;
      s = h/z;
      f = x*c + g*s;
      g = - x*s + g*c;
      h = y*s;
      y *= c;
      for (int j=0; j<n; j++) {
        x = ab[j+(i-1)*n];
        z = ab[j+i*n];
        ab[j+(i-1)*n] = x*c + z*s;
        ab[j+i*n] = - x*s + z*c;
      }
      if (std::fabs(f) < std::fabs(h)) {
	      FloatType fh = f/h;
	      z = std::fabs(h)*std::sqrt(1.0+fh*fh);
      } else {
	      FloatType hf = h/f;
	      z = (f!=0.0 ? std::fabs(f)*std::sqrt(1.0+hf*hf) : 0.0);
      }
      q[i-1] = z;
      if (z == 0.0) z = f = 1.0;
      c = f/z;
      s = h/z;
      f = c*g + s*y;
      x = - s*g + c*y;
    }
    e[l] = 0.0;
    e[k] = f;
    q[k] = x;
    goto TestFsplitting;
Convergence:
    if (z < 0.0){
      q[k] = - z;
      for (int j=0; j<n; j++) ab[j+k*n] = - ab[j+k*n];
    }
  }
}

} // namespace PRAXIS
