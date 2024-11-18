// Pure Go implementation of the Matrix arithmetic and GSL libs
// draft version
// unstable API
//
// Statistics(GSL)
//
// Adjusted for concurrency and matrices specifics GSL Statistics library
// for simple vectors you should use GSL lib because it's more camfotable
// [github.com/grd/stat]
// stat lib
// new: https://github.com/grd/stat
// old: https://github.com/grd/statistics
package matrix_arithmetic

// TODO: Unified Interface for v.2

import (
	// common purpose
	_ "bufio"
	_ "errors" // errors.New()
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

// Mean calculates the arithmetic mean with the recurrence relation
func Mean(a *Matrix, cursor int, vector uint8) (mean float64) {
	var n int
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			mean += (value - mean) / float64(j+1)
			//		fmt.Printf("[%d][%d] value: %.2f     mean:%.2f\n", cursor, j, value, mean)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			mean += (value - mean) / float64(i+1)
			//		fmt.Printf("[%d][%d] value: %.2f     mean:%.2f\n", cursor, i, value, mean)
		}
	}

	return
}

// Mean square calculates the arithmetic mean of squeare (x*x) with the recurrence relation
func MeanSquare(a, b *Matrix, cursor_a, cursor_b int, vector uint8) (meanSquare float64) {
	var n int
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < a.rows; j++ {
			value_a, _ := a.Getij(cursor_a, j)
			value_b, _ := b.Getij(cursor_b, j)
			meanSquare += (value_a*value_b - meanSquare) / float64(j+1)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_a, _ := a.Getij(i, cursor_a)
			value_b, _ := b.Getij(i, cursor_b)
			meanSquare += (value_a*value_b - meanSquare) / float64(i+1)
		}
	}

	return
}

func Absdev(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	return AbsdevMean(a, cursor, vector, mean)
}

// AbsdevMean finds the absolute deviation of the data interface
func AbsdevMean(a *Matrix, cursor int, vector uint8, mean float64) float64 {
	var sum float64 = float64(0)
	var n int = 0

	// the sum of the absolute deviations
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			sum += math.Abs(value - mean)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			sum += math.Abs(value - mean)
		}
	}

	return sum / float64(n)
}

// ver form DSL
// takes a dataset and calculates the covariance
/*
func covariance(a, b *Matrix, cursor_a, cursor_b int, vector uint8, mean_a, mean_b float64) (covar float64) {
	var n int
	// calculate the sum of the squares
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value_a,_ := a.Getij(cursor_a,j)
			value_b,_ := b.Getij(cursor_b,j)
			delta_a := (value_a - mean_a)
			delta_b := (value_b - mean_b)
			covar += (delta_a * delta_b - covar) / float64(j+1)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_a,_ := a.Getij(i,cursor_a)
			value_b,_ := b.Getij(i,cursor_b)
			delta_a := (value_a - mean_a)
			delta_b := (value_b - mean_b)
			covar += (delta_a * delta_b - covar) / float64(i+1)
		}
	}

	return
}
*/

// my ver
func covariance(a, b *Matrix, cursor_a, cursor_b int, vector uint8, mean_a, mean_b float64) (covar float64) {
	meanSquare := MeanSquare(a, b, cursor_a, cursor_b, vector)
	// covarince as COV[xy] = M[x*y] - M[x]*M[y]
	covar = meanSquare - mean_a*mean_b
	return
}

// ver form DSL
/*
func CovarianceMean(a, b *Matrix, cursor_a, cursor_b int, vector uint8, mean_a, mean_b float64) float64 {
	var n int
	var covar float64

	if vector == ROWVECTOR {
		n = a.cols
		covar = covariance(a, b, cursor_a, cursor_b, vector, mean_a, mean_b)
	} else {
		n = a.rows
		covar = covariance(a, b, cursor_a, cursor_b, vector, mean_a, mean_b)
	}

	return covar * float64(n) / float64(n-1)
}
*/

// my ver
func CovarianceMean(a, b *Matrix, cursor_a, cursor_b int, vector uint8, mean_a, mean_b float64) float64 {
	return covariance(a, b, cursor_a, cursor_b, vector, mean_a, mean_b)
}

func Covariance(a, b *Matrix, cursor_a, cursor_b int, vector uint8) float64 {
	mean_a := Mean(a, cursor_a, vector)
	mean_b := Mean(b, cursor_b, vector)

	return CovarianceMean(a, b, cursor_a, cursor_b, vector, mean_a, mean_b)
}

// Correlation()
//
//	Calculate Pearson correlation = cov(X, Y) / (sigma_X * sigma_Y)
//
// This routine efficiently computes the correlation in one pass of the
// data and makes use of the algorithm described in:
// B. P. Welford, "Note on a Method for Calculating Corrected Sums of
// Squares and Products", Technometrics, Vol 4, No 3, 1962.
// This paper derives a numerically stable recurrence to compute a sum
// of products
// S = sum_{i=1..N} [ (x_i - mu_x) * (y_i - mu_y) ]
// with the relation
// S_n = S_{n-1} + ((n-1)/n) * (x_n - mu_x_{n-1}) * (y_n - mu_y_{n-1})
func Correlation(a, b *Matrix, cursor_a, cursor_b int, vector uint8) (correl float64) {
	var sum_asq, sum_bsq, sum_cross float64
	var n int

	//
	// Compute:
	// sum_xsq = Sum [ (x_i - mu_x)^2 ],
	// sum_ysq = Sum [ (y_i - mu_y)^2 ] and
	// sum_cross = Sum [ (x_i - mu_x) * (y_i - mu_y) ]
	// using the above relation from Welford's paper
	//

	if vector == ROWVECTOR {
		n = a.cols

		value_a, _ := a.Getij(cursor_a, 0)
		value_b, _ := b.Getij(cursor_b, 0)
		mean_a := value_a
		mean_b := value_b

		for j := 1; j < n; j++ {
			ratio := float64(j) / float64(j+1)
			value_a, _ = a.Getij(cursor_a, j)
			value_b, _ = b.Getij(cursor_b, j)
			delta_a := value_a - mean_a
			delta_b := value_b - mean_b
			sum_asq += delta_a * delta_a * ratio
			sum_bsq += delta_b * delta_b * ratio
			sum_cross += delta_a * delta_b * ratio
			mean_a += delta_a / float64(j+1)
			mean_b += delta_b / float64(j+1)
		}

	} else {
		n = a.rows

		value_a, _ := a.Getij(0, cursor_a)
		value_b, _ := b.Getij(0, cursor_b)
		mean_a := value_a
		mean_b := value_b

		for i := 1; i < n; i++ {
			ratio := float64(i) / float64(i+1)
			value_a, _ = a.Getij(i, cursor_a)
			value_b, _ = b.Getij(i, cursor_b)
			delta_a := value_a - mean_a
			delta_b := value_b - mean_b
			sum_asq += delta_a * delta_a * ratio
			sum_bsq += delta_b * delta_b * ratio
			sum_cross += delta_a * delta_b * ratio
			mean_a += delta_a / float64(i+1)
			mean_b += delta_b / float64(i+1)
		}

	}

	correl = sum_cross / (math.Sqrt(sum_asq) * math.Sqrt(sum_bsq))
	return
}

// ver from DSL
/*
func _variance(a *Matrix, cursor int, vector uint8, mean float64) (variance float64) {
	var n int

	// calculate the sum of the squares
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value,_ := a.Getij(cursor,j)
			delta := value - mean
			// TODO: long double for variance... How to implement in Go?
			variance += ((delta * delta) - variance) / float64(j+1)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value,_ := a.Getij(i,cursor)
			delta := value - mean
			// TODO: long double for variance... How to implement in Go?
			variance += ((delta * delta) - variance) / float64(i+1)
		}
	}

	return
}
*/

// my ver
func _variance(a *Matrix, cursor int, vector uint8, mean float64) (variance float64) {
	meanSquare := MeanSquare(a, a, cursor, cursor, vector)
	//variance as D[x] = M[x^2] - M^2[x]
	return meanSquare - mean*mean
}

func VarianceWithFixedMean(a *Matrix, cursor int, vector uint8, mean float64) float64 {
	return _variance(a, cursor, vector, mean)
}

func SdWithFixedMean(a *Matrix, cursor int, vector uint8, mean float64) float64 {
	variance := _variance(a, cursor, vector, mean)
	return math.Sqrt(variance)
}

func VarianceMean(a *Matrix, cursor int, vector uint8, mean float64) float64 {
	var n int
	variance := _variance(a, cursor, vector, mean)

	if vector == ROWVECTOR {
		n = a.cols
	} else {
		n = a.rows
	}

	return variance * float64(n) / float64(n-1)
}

func SdMean(a *Matrix, cursor int, vector uint8, mean float64) float64 {
	var n int

	variance := _variance(a, cursor, vector, mean)

	if vector == ROWVECTOR {
		n = a.cols
	} else {
		n = a.rows
	}

	return math.Sqrt(variance * float64(n) / float64(n-1))
}

// ver from DSL
//func Variance(a *Matrix, cursor int, vector uint8) float64 {
//	mean := Mean(a, cursor, vector)
//	return VarianceMean(a, cursor, vector, mean)
//}

// my versions
// VAR.P in calc
// Calculates a variance based on the entire population.
func VarianceP(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	return CovarianceMean(a, a, cursor, cursor, vector, mean, mean)
}

// VAR in calc
// Estimates the variance based on a sample.
func Variance(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	return VarianceMean(a, cursor, vector, mean)
}

func Sd(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	return SdMean(a, cursor, vector, mean)
}

// TssMean takes a dataset and finds the sum of squares about the mean
func TssMean(a *Matrix, cursor int, vector uint8, mean float64) (res float64) {
	var n int
	// find the sum of the squares
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			delta := value - mean
			res += delta * delta
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			delta := value - mean
			res += delta * delta
		}
	}

	return
}

func Tss(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	return TssMean(a, cursor, vector, mean)
}

func Kurtosis(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	est_sd := SdMean(a, cursor, vector, mean)
	return KurtosisMainSd(a, cursor, vector, mean, est_sd)
}

func KurtosisMainSd(a *Matrix, cursor int, vector uint8, mean, sd float64) float64 {
	var avg, kurtosis float64
	var n int

	// calculate the fourth moment the deviations, normalized by the sd

	// we use a recurrence relation to stable update a running value so
	// there aren't any large sums that can overflow

	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			x := (value - mean) / sd
			avg += (x*x*x*x - avg) / float64(j+1)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			x := (value - mean) / sd
			avg += (x*x*x*x - avg) / float64(i+1)
		}
	}

	kurtosis = avg - 3.0 // makes kurtosis zero for a Gaussian

	return kurtosis
}

func Lag1Autocorrelation(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	return Lag1AutocorrelationMean(a, cursor, vector, mean)
}

func Lag1AutocorrelationMean(a *Matrix, cursor int, vector uint8, mean float64) float64 {
	var r1, q, v float64
	var n int

	if vector == ROWVECTOR {
		n = a.cols

		value, _ := a.Getij(cursor, 0)
		v = (value - mean) * (value - mean)

		for j := 1; j < n; j++ {
			value0, _ := a.Getij(cursor, j-1)
			value1, _ := a.Getij(cursor, j)
			delta0 := value0 - mean
			delta1 := value1 - mean
			q += (delta0*delta1 - q) / float64(j+1)
			v += (delta1*delta1 - v) / float64(j+1)
		}

	} else {
		n = a.rows

		value, _ := a.Getij(0, cursor)
		v = (value - mean) * (value - mean)

		for i := 1; i < n; i++ {
			value0, _ := a.Getij(i-1, cursor)
			value1, _ := a.Getij(i, cursor)
			delta0 := value0 - mean
			delta1 := value1 - mean
			q += (delta0*delta1 - q) / float64(i+1)
			v += (delta1*delta1 - v) / float64(i+1)
		}
	}

	r1 = q / v

	return r1
}

// MedianFromSortedData calculates the median of the sorted data.
// Note that the function doesn't check wheather the data is actually sorted.
func MedianFromSortedData(a *Matrix, cursor int, vector uint8) (median float64) {
	var n int

	if vector == ROWVECTOR {
		n = a.cols

		if n == 0 {
			return float64(0)
		}

		lhs := (n - 1) / 2
		rhs := n / 2

		if lhs == rhs {
			value, _ := a.Getij(cursor, lhs)
			median = value
		} else {
			value_lhs, _ := a.Getij(cursor, lhs)
			value_rhs, _ := a.Getij(cursor, rhs)
			median = (value_lhs + value_rhs) / float64(2)
		}

	} else {
		n = a.rows

		if n == 0 {
			return float64(0)
		}

		lhs := (n - 1) / 2
		rhs := n / 2

		if lhs == rhs {
			value, _ := a.Getij(cursor, lhs)
			median = value
		} else {
			value_lhs, _ := a.Getij(lhs, cursor)
			value_rhs, _ := a.Getij(rhs, cursor)
			median = (value_lhs + value_rhs) / float64(2)
		}
	}

	return
}

// Max finds the first largest member and the members position within the data
func Max(a *Matrix, cursor int, vector uint8) (max float64, max_index int) {
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		value, _ := a.Getij(cursor, 0)
		max = value

		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			xj := value

			if xj > max {
				max = xj
				max_index = j
			}

			if math.IsNaN(xj) {
				max = xj
				max_index = j
				return
			}
		}

	} else {
		n = a.rows
		value, _ := a.Getij(0, cursor)
		max = value

		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			xi := value

			if xi > max {
				max = xi
				max_index = i
			}

			if math.IsNaN(xi) {
				max = xi
				max_index = i
				return
			}
		}

	}

	return
}

// Min finds the first smallest member and the members position within the data
func Min(a *Matrix, cursor int, vector uint8) (min float64, min_index int) {
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		value, _ := a.Getij(cursor, 0)
		min = value

		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			xj := value

			if xj < min {
				min = xj
				min_index = j
			}

			if math.IsNaN(xj) {
				min = xj
				min_index = j
				return
			}
		}

	} else {
		n = a.rows
		value, _ := a.Getij(0, cursor)
		min = value

		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			xi := value

			if xi < min {
				min = xi
				min_index = i
			}

			if math.IsNaN(xi) {
				min = xi
				min_index = i
				return
			}
		}

	}

	return
}

// Minmax finds the first smallest and largest members and
// the members positions within the data
func Minmax(a *Matrix, cursor int, vector uint8) (min float64, min_index int, max float64, max_index int) {
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		value, _ := a.Getij(cursor, 0)
		max = value
		min = value

		stop := false
		for j := 0; j < n && !stop; j++ {
			value, _ := a.Getij(cursor, j)
			xj := value

			if xj < min {
				min = xj
				min_index = j
			}

			if xj > max {
				max = xj
				max_index = j
			}

			if math.IsNaN(xj) {
				min = xj
				max = xj
				min_index = j
				max_index = j
				stop = true
			}
		}

	} else {
		n = a.rows
		value, _ := a.Getij(0, cursor)
		max = value
		min = value

		stop := false
		for i := 0; i < n && !stop; i++ {
			value, _ := a.Getij(i, cursor)
			xi := value

			if xi < min {
				min = xi
				min_index = i
			}

			if xi > max {
				max = xi
				max_index = i
			}

			if math.IsNaN(xi) {
				min = xi
				max = xi
				min_index = i
				max_index = i
				stop = true
			}
		}

	}

	return
}

// PVariance finds the pooled variance of two datasets
func PVariance(a, b *Matrix, cursor_a, cursor_b int, vector uint8) float64 {
	var n_a int
	var n_b int
	var var_a float64
	var var_b float64

	if vector == ROWVECTOR {
		n_a = a.cols
		n_b = b.cols
	} else {
		n_a = a.rows
		n_b = b.rows
	}

	var_a = Variance(a, cursor_a, vector)
	var_b = Variance(b, cursor_b, vector)

	return ((float64(n_a-1) * var_a) + (float64(n_b-1) * var_b)) / float64(n_a+n_b-2)
}

// QuantileFromSortedData performs the quantile function, also called percent
// point function or inverse cumulative distribution function, on the sorted data.
// Note that the function doesn't check wheather the data is actually sorted.
func QuantileFromSortedData(a *Matrix, cursor int, vector uint8, f float64) (result float64) {

	var n int

	if vector == ROWVECTOR {
		n = a.cols
		index := f * float64(n-1)
		lhs := int(index)
		delta := index - float64(lhs)

		if lhs == n-1 {
			value, _ := a.Getij(cursor, lhs)
			result = value
		} else {
			value1, _ := a.Getij(cursor, lhs)
			value2, _ := a.Getij(cursor, lhs+1)
			result = (1-delta)*value1 + delta*value2
		}
	} else {
		n = a.rows
		index := f * float64(n-1)
		lhs := int(index)
		delta := index - float64(lhs)

		if lhs == n-1 {
			value, _ := a.Getij(lhs, cursor)
			result = value
		} else {
			value1, _ := a.Getij(lhs, cursor)
			value2, _ := a.Getij(lhs+1, cursor)
			result = (1-delta)*value1 + delta*value2
		}
	}

	return
}

func Skew(a *Matrix, cursor int, vector uint8) float64 {
	mean := Mean(a, cursor, vector)
	sd := SdMean(a, cursor, vector, mean)
	return SkewMeanSd(a, cursor, vector, mean, sd)
}

// SkewMeanSd calculates the skewness of a dataset
func SkewMeanSd(a *Matrix, cursor int, vector uint8, mean, sd float64) (skew float64) {
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			x := (value - mean) / sd
			skew += (x*x*x - skew) / float64(j+1)
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			x := (value - mean) / sd
			skew += (x*x*x - skew) / float64(i+1)
		}
	}

	return
}

func WAbsdev(w, a *Matrix, cursor_w, cursor_a int, vector uint8) float64 {
	wmean := WMean(w, a, cursor_w, cursor_a, vector)
	return WAbsdevMean(w, a, cursor_w, cursor_a, vector, wmean)
}

// WAbsdevMean calculates the weighted absolute deviation of a dataset
func WAbsdevMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) (wabsdev float64) {
	var W float64
	var n int

	// calculate the sum of the absolute deviations
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value_w, _ := w.Getij(cursor_w, j)
			wj := value_w

			if wj > 0 {
				value_a, _ := a.Getij(cursor_a, j)
				delta := math.Abs(value_a - wmean)
				W += wj
				wabsdev += (delta - wabsdev) * (wj / W)
			}
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_w, _ := w.Getij(i, cursor_w)
			wi := value_w

			if wi > 0 {
				value_a, _ := a.Getij(i, cursor_a)
				delta := math.Abs(value_a - wmean)
				W += wi
				wabsdev += (delta - wabsdev) * (wi / W)
			}
		}
	}

	return
}

func WKurtosis(w, a *Matrix, cursor_w, cursor_a int, vector uint8) float64 {
	wmean := WMean(w, a, cursor_w, cursor_a, vector)
	wsd := WSdMean(w, a, cursor_w, cursor_a, vector, wmean)
	return WKurtosisMeanSd(w, a, cursor_w, cursor_a, vector, wmean, wsd)
}

// WKurtosisMean calculates the kurtosis of a dataset
func WKurtosisMeanSd(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean, wsd float64) float64 {
	var wavg, W float64
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value_w, _ := w.Getij(cursor_w, j)
			wj := value_w

			if wj > 0 {
				value_a, _ := a.Getij(cursor_a, j)
				x := (value_a - wmean) / wsd
				W += wj
				wavg += (x*x*x*x - wavg) * (wj / W)
			}
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_w, _ := w.Getij(i, cursor_w)
			wi := value_w

			if wi > 0 {
				value_a, _ := a.Getij(i, cursor_a)
				x := (value_a - wmean) / wsd
				W += wi
				wavg += (x*x*x*x - wavg) * (wi / W)
			}
		}
	}

	return wavg - 3.0 // makes kurtosis zero for a Gaussian
}

func wvariance(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) (wvariance float64) {
	var W float64
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value_w, _ := w.Getij(cursor_w, j)
			wj := value_w

			if wj > 0 {
				value_a, _ := a.Getij(cursor_a, j)
				delta := value_a - wmean
				W += wj
				wvariance += (delta*delta - wvariance) * (wj / W)
			}
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_w, _ := w.Getij(i, cursor_w)
			wi := value_w

			if wi > 0 {
				value_a, _ := a.Getij(i, cursor_a)
				delta := value_a - wmean
				W += wi
				wvariance += (delta*delta - wvariance) * (wi / W)
			}
		}
	}

	return
}

func factor(w *Matrix, cursor int, vector uint8) (factor float64) {
	var a, b float64
	var n int

	// the sum of the squares
	if vector == ROWVECTOR {
		n = w.cols
		for j := 0; j < n; j++ {
			value, _ := w.Getij(cursor, j)
			wj := value

			if wj > 0 {
				a += wj
				b += wj * wj
			}
		}
	} else {
		n = w.rows
		for i := 0; i < n; i++ {
			value, _ := w.Getij(i, cursor)
			wi := value

			if wi > 0 {
				a += wi
				b += wi * wi
			}
		}
	}

	factor = (a * a) / ((a * a) - b)

	return
}

func WVarianceWithFixedMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) float64 {
	return wvariance(w, a, cursor_w, cursor_a, vector, wmean)
}

func WsdWithFixedMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) float64 {
	wvariance := wvariance(w, a, cursor_w, cursor_a, vector, wmean)
	return math.Sqrt(wvariance)
}

func WVarianceMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) float64 {
	variance := wvariance(w, a, cursor_w, cursor_a, vector, wmean)
	scale := factor(w, cursor_w, vector)

	return scale * variance
}

func WSdMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) float64 {
	variance := wvariance(w, a, cursor_w, cursor_a, vector, wmean)
	scale := factor(w, cursor_w, vector)
	return math.Sqrt(scale * variance)
}

func WSd(w, a *Matrix, cursor_w, cursor_a int, vector uint8) float64 {
	wmean := WMean(w, a, cursor_w, cursor_a, vector)
	return WSdMean(w, a, cursor_w, cursor_a, vector, wmean)
}

func WVariance(w, a *Matrix, cursor_w, cursor_a int, vector uint8) float64 {
	wmean := WMean(w, a, cursor_w, cursor_a, vector)
	return WVarianceMean(w, a, cursor_w, cursor_a, vector, wmean)
}

// WTssMean takes a dataset and finds the weighted sum of squares about wmean
func WTssMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8, wmean float64) (res float64) {
	var n int

	// find the sum of the squares
	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value_w, _ := w.Getij(cursor_w, j)
			wj := value_w

			if wj > 0 {
				value_a, _ := a.Getij(cursor_a, j)
				delta := value_a - wmean
				res += wj * delta * delta
			}
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_w, _ := w.Getij(i, cursor_w)
			wi := value_w

			if wi > 0 {
				value_a, _ := a.Getij(i, cursor_a)
				delta := value_a - wmean
				res += wi * delta * delta
			}
		}
	}

	return
}

func WTss(w, a *Matrix, cursor_w, cursor_a int, vector uint8) float64 {
	wmean := WMean(w, a, cursor_w, cursor_a, vector)
	return WTssMean(w, a, cursor_w, cursor_a, vector, wmean)
}

// WMean calculates the weighted arithmetic mean of a dataset
func WMean(w, a *Matrix, cursor_w, cursor_a int, vector uint8) (wmean float64) {
	var W float64
	var n int

	if vector == ROWVECTOR {
		n = a.cols
		for j := 0; j < n; j++ {
			value_w, _ := w.Getij(cursor_w, j)
			wj := value_w

			if wj > 0 {
				W += wj
				value_a, _ := a.Getij(cursor_a, j)
				wmean += (value_a - wmean) * (wj / W)
			}
		}
	} else {
		n = a.rows
		for i := 0; i < n; i++ {
			value_w, _ := w.Getij(i, cursor_w)
			wi := value_w

			if wi > 0 {
				W += wi
				value_a, _ := a.Getij(i, cursor_a)
				wmean += (value_a - wmean) * (wi / W)
			}
		}
	}

	return
}
