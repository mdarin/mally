// Pure Go implementation of the singular-value decomposition (SVD)
// draft version
// unstable API
//
// # Linear algebra
//
// SVD is a factorization of a real or complex matrix.
// It is the generalization of the eigendecomposition of a positive semidefinite
// normal matrix (for example, a symmetric matrix with positive eigenvalues)
// to any m × n matrix via an extension of the polar decomposition.
// It has many useful applications in signal processing and statistics.
package matrix_arithmetic

import (
	// common purpose
	_ "bufio"
	"errors" // errors.New()
	_ "fmt"
	_ "io"
	_ "log"
	"math"
	_ "os"
	_ "strconv" // aoti and so on convertions
	_ "strings"
	_ "time"
	// "sync"
)

const ()

var ()

/*
 * SVD - singular-value decomposition routine.
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
 * code from Numerical Recipes adapted by Luke Tierney and David Betz
 * translated to Go by Michael Darin
 *
 * Input to SVD is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a [internal]
 *   n = column dimension of a [internal]
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
 * Output
 *   err = error desc or nil on success
 */
func SVD(a, w, v *Matrix) error {
	// declare variables
	var flag int
	var i int
	var its int
	var j int
	var jj int
	var k int
	var l int
	var nm int
	var c float64
	var f float64
	var h float64
	var s float64
	var x float64
	var y float64
	var z float64
	var anorm float64 = 0.0
	var g float64 = 0.0
	var scale float64 = 0.0

	m := a.Rows()
	n := a.Cols()

	if m < n {
		return errors.New("#rows must be more or equal #cols")
	}

	// vector
	rv1, _ := New(n, 1)

	// Householder reduction to bidiagonal form
	for i = 0; i < n; i++ {

		// left-hand reduction
		l = i + 1
		rv1.Setij(i, 0, scale*g)
		g = 0.0
		s = 0.0
		scale = 0.0
		if i < m {
			for k = i; k < m; k++ {
				value_a, _ := a.Getij(k, i)
				scale += math.Abs(value_a)
			}
			if scale != 0.0 {
				for k = i; k < m; k++ {
					value_a, _ := a.Getij(k, i)
					ratio := value_a / scale
					a.Setij(k, i, ratio)
					s += ratio * ratio
				}
				f, _ = a.Getij(i, i)
				g = -sign(math.Sqrt(s), f)
				h = f*g - s
				a.Setij(i, i, f-g)
				if i != n-1 {
					for j = l; j < n; j++ {
						s = 0.0
						for k = i; k < m; k++ {
							value_ai, _ := a.Getij(k, i)
							value_aj, _ := a.Getij(k, j)
							s += value_ai * value_aj
						}
						f = s / h
						for k = i; k < m; k++ {
							value_ai, _ := a.Getij(k, i)
							value_aj, _ := a.Getij(k, j)
							value_aj += f * value_ai
							a.Setij(k, j, value_aj)
						}
					}
				} // eof if i != n-1
				for k = i; k < m; k++ {
					value_a, _ := a.Getij(k, i)
					a.Setij(k, i, value_a*scale)
				}
			} //eof if scale
		} //eof if i<m
		w.Setij(i, 0, scale*g)

		// right-hand reduction
		g = 0.0
		s = 0.0
		scale = 0.0
		if i < m && i != n-1 {
			for k = l; k < n; k++ {
				value_a, _ := a.Getij(i, k)
				scale += math.Abs(value_a)
			}
			if scale != 0.0 {
				for k = l; k < n; k++ {
					value_a, _ := a.Getij(i, k)
					ratio := value_a / scale
					a.Setij(i, k, ratio)
					s += ratio * ratio
				}
				f, _ = a.Getij(i, l)
				g = -sign(math.Sqrt(s), f)
				h = f*g - s
				a.Setij(i, l, f-g)
				for k = l; k < n; k++ {
					value_a, _ := a.Getij(i, k)
					rv1.Setij(k, 0, value_a/h)
				}
				if i != m-1 {
					for j = l; j < m; j++ {
						s = 0.0
						for k = l; k < n; k++ {
							value_aj, _ := a.Getij(j, k)
							value_ai, _ := a.Getij(i, k)
							s += value_aj * value_ai
						}
						for k = l; k < n; k++ {
							value_a, _ := a.Getij(j, k)
							value_rv1, _ := rv1.Getij(k, 0)
							value_a += s * value_rv1
							a.Setij(j, k, value_a)
						}
					} // eof for j=l
				} // eof if i!=m-1
				for k = l; k < n; k++ {
					value_a, _ := a.Getij(i, k)
					a.Setij(i, k, value_a*scale)
				}
			} // eof if scale
		} // eof if i<m i!=n-1
		value_w, _ := w.Getij(i, 0)
		value_rv1, _ := rv1.Getij(i, 0)
		anorm = math.Max(anorm, (math.Abs(value_w) + math.Abs(value_rv1)))
	}

	// accumulate the right-hand transformation
	for i = n - 1; i >= 0; i-- {
		if i < n-1 {
			if g != 0.0 {
				// double division to avoid underflow
				for j = l; j < n; j++ {
					value_aj, _ := a.Getij(i, j)
					value_al, _ := a.Getij(i, l)
					v.Setij(j, i, value_aj/value_al/g)
				}
				for j = l; j < n; j++ {
					s = 0.0
					for k = l; k < n; k++ {
						value_a, _ := a.Getij(i, k)
						value_v, _ := v.Getij(k, j)
						s += value_a * value_v
					}
					for k = l; k < n; k++ {
						value_vj, _ := v.Getij(k, j)
						value_vi, _ := v.Getij(k, i)
						value_vj += s * value_vi
						v.Setij(k, j, value_vj)
					}
				} // eof for j=l
			} // eof g!=0
			for j = l; j < n; j++ {
				v.Setij(i, j, 0.0)
				v.Setij(j, i, 0.0)
			}
		} // eof if i<n-1
		v.Setij(i, i, 1.0)
		g, _ = rv1.Getij(i, 0)
		l = i
	} // eof for i=n-1

	// accumulate the left-hand transformation
	for i = n - 1; i >= 0; i-- {
		l = i + 1
		g, _ = w.Getij(i, 0)
		if i < n-1 {
			for j = l; j < n; j++ {
				a.Setij(i, j, 0.0)
			}
		}
		if g != 0.0 {
			g = 1.0 / g
			if i != n-1 {
				for j = l; j < n; j++ {
					s = 0.0
					for k = l; k < m; k++ {
						value_ai, _ := a.Getij(k, i)
						value_aj, _ := a.Getij(k, j)
						s += value_ai * value_aj
					}
					value_a, _ := a.Getij(i, i)
					f = s / value_a * g
					for k = i; k < m; k++ {
						value_aj, _ := a.Getij(k, j)
						value_ai, _ := a.Getij(k, i)
						value_aj += f * value_ai
						a.Setij(k, j, value_aj)
					}
				} // eof for j=l
			} // eof if i!=n-1
			for j = i; j < m; j++ {
				value_a, _ := a.Getij(j, i)
				a.Setij(j, i, value_a*g)
			}
		} else {
			for j = i; j < m; j++ {
				a.Setij(j, i, 0.0)
			}
		} // eof if g!=0
		value_a, _ := a.Getij(i, i)
		value_a++
		a.Setij(i, i, value_a)
	} // eof for i=n-1

	// diagonalize the bidiagonal form
	for k = n - 1; k >= 0; k-- { // loop over singular values
		for its = 0; its < 30; its++ { // loop over allowed iterations
			flag = 1
			for l = k; l >= 0; l-- { // test for splitting
				nm = l - 1
				value_rv1, _ := rv1.Getij(l, 0)
				if math.Abs(value_rv1)+anorm == anorm {
					flag = 0
					break
				}
				value_w, _ := w.Getij(nm, 0)
				if math.Abs(value_w)+anorm == anorm {
					break
				}
			} // eof for l=k
			if flag != 0 {
				c = 0.0
				s = 1.0
				for i = l; i <= k; i++ {
					value_rv1, _ := rv1.Getij(i, 0)
					f = s * value_rv1
					if math.Abs(f)+anorm != anorm {
						g, _ = w.Getij(i, 0)
						h = pythag(f, g)
						w.Setij(i, 0, h)
						h = 1.0 / h
						c = g * h
						s = (-f * h)
						for j = 0; j < m; j++ {
							y, _ = a.Getij(j, nm)
							z, _ = a.Getij(j, i)
							a.Setij(j, nm, (y*c + z*s))
							a.Setij(j, i, (z*c - y*s))
						}
					} // eof math.Abs()
				} // eof for i=l
			} // eof if flag
			z, _ = w.Getij(k, 0)
			if l == k { // convergence
				if z < 0.0 { // make singular value nonnegative
					w.Setij(k, 0, (-z))
					for j = 0; j < n; j++ {
						value_v, _ := v.Getij(j, k)
						v.Setij(j, k, (-value_v))
					}
				}
				break
			}
			if its >= 30 {
				return errors.New("No convergence after 30 iterations")
			}

			// shift from bottom 2 x 2 minor
			x, _ = w.Getij(l, 0)
			nm = k - 1
			y, _ = w.Getij(nm, 0)
			g, _ = rv1.Getij(nm, 0)
			h, _ = rv1.Getij(k, 0)
			f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0 * h * y)
			g = pythag(f, 1.0)
			f = ((x-z)*(x+z) + h*((y/(f+sign(g, f)))-h)) / x

			// next QR transformation
			c = 1.0
			s = 1.0
			for j = l; j <= nm; j++ {
				i = j + 1
				g, _ = rv1.Getij(i, 0)
				y, _ = w.Getij(i, 0)
				h = s * g
				g = c * g
				z = pythag(f, h)
				rv1.Setij(i, 0, z)
				c = f / z
				s = h / z
				f = x*c + g*s
				g = g*c - x*s
				h = y * s
				y = y * c
				for jj = 0; jj < n; jj++ {
					x, _ = v.Getij(jj, j)
					z, _ = v.Getij(jj, i)
					v.Setij(jj, j, (x*c + z*s))
					v.Setij(jj, i, (z*c - x*s))
				}
				z = pythag(f, h)
				w.Setij(j, 0, z)
				if z != 0.0 {
					z = 1.0 / z
					c = f * z
					s = h * z
				}
				f = (c * g) + (s * y)
				x = (c * y) - (s * g)
				for jj = 0; jj < m; jj++ {
					y, _ = a.Getij(jj, j)
					z, _ = a.Getij(jj, i)
					a.Setij(jj, j, (y*c + z*s))
					a.Setij(jj, i, (z*c - y*s))
				}
			} // for j=l
			rv1.Setij(l, 0, 0.0)
			rv1.Setij(k, 0, f)
			w.Setij(k, 0, x)
		} // eof loop over allowed iterations
	} // eof loop over singular values
	return nil
} // eof SVD procedure

/**
 ** Internals
 */

/* Computes (a^2 + b^2 )^1/2 without destructive underflow or overflow. */
func pythag(a, b float64) float64 {
	at := math.Abs(a)
	bt := math.Abs(b)
	ct := 0.0
	result := 0.0

	if at > bt {
		ct = bt / at
		result = at * math.Sqrt(1.0+ct*ct)
	} else if bt > 0.0 {
		ct = at / bt
		result = bt * math.Sqrt(1.0+ct*ct)
	} //else result = 0.0

	return result
}

/* The word sign refers to the property of being positive or negative. */
func sign(a, b float64) float64 {
	if b >= 0.0 {
		return math.Abs(a)
	} else {
		return -math.Abs(a)
	}
}
