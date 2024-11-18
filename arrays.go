// Pure Go implementation of the Matrix arithmetic and GSL libs
// draft version
// unstable API
//
// # Array Functions and Mathematical Functions
//
// Calc functions involved in interaction with this matrix arithmetic library
// This module provides a set of handy array and math functions
package matrix_arithmetic

// FIXME: EXPERIMENTAL!
// TODO: Unified Interface for v.2

import (
	// common purpose
	_ "bufio"
	_ "errors" // errors.New()
	_ "fmt"
	_ "io"
	_ "log"
	_ "math"
	_ "os"
	_ "strconv" // aoti and so on convertions
	_ "strings"
	_ "time"

	// simple concurrency
	_ "sync"
)

/*
SUBTOTAL
Calculates subtotals. If a range already contains subtotals, these are not used for further calculations. Use this function with the AutoFilters to take only the filtered records into account.

Syntax
SUBTOTAL(Function; Range)

Function is a number that stands for one of the following functions:
Function index / Function
1 AVERAGE
2 COUNT
3 COUNTA
4 MAX
5 MIN
6 PRODUCT
7 STDEV
8 STDEVP
9 SUM
10 VAR
11 VARP

Range is the range whose cells are included.

Example
You have a table in the cell range A1:B5 containing cities in column A and accompanying figures in column B. You have used an AutoFilter so that you only see rows containing the city Hamburg. You want to see the sum of the figures that are displayed; that is, just the subtotal for the filtered rows. In this case the correct formula would be:
=SUBTOTAL(9;B2:B5)
*/

/*
SUM
Adds all the numbers in a range of cells.

Syntax
SUM(Number1; Number2; ...; Number30)
Number 1 to Number 30 are up to 30 arguments whose sum is to be calculated.

Example
If you enter the numbers 2; 3 and 4 in the Number 1; 2 and 3 text boxes, 9 will be returned as the result.
=SUM(A1;A3;B5) calculates the sum of the three cells. =SUM (A1:E10) calculates the sum of all cells in the A1 to E10 cell range.
Conditions linked by AND can be used with the function SUM() in the following manner:
Example assumption: You have entered invoices into a table. Column A contains the date value of the invoice, column B the amounts. You want to find a formula that you can use to return the total of all amounts only for a specific month, e.g. only the amount for the period >=2008-01-01 to <2008-02-01. The range with the date values covers A1:A40, the range containing the amounts to be totaled is B1:B40. C1 contains the start date, 2008-01-01, of the invoices to be included and C2 the date, 2008-02-01, that is no longer included.
Enter the following formula as an array formula:
=SUM((A1:A40>=C1)*(A1:A40<C2)*B1:B40)
In order to enter this as an array formula, you must press the Shift+ Ctrl+ Enter keys instead of simply pressing the Enter key to close the formula. The formula will then be shown in the Formula bar enclosed in braces.
{=SUM((A1:A40>=C1)*(A1:A40<C2)*B1:B40)}
The formula is based on the fact that the result of a comparison is 1 if the criterion is met and 0 if it is not met. The individual comparison results will be treated as an array and used in matrix multiplication, and at the end the individual values will be totaled to give the result matrix.
*/

// sum of all table cells
func SumM(a *Matrix) (sum float64) {
	return PsumM(a, 0, 0, 0, 0)
}

// partial sum table
func PsumM(a *Matrix, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	//index := lo_i * lo_j
	for i := 0 + lo_i; i < a.Rows()-hi_i; i++ {
		for j := 0 + lo_j; j < a.Cols()-hi_j; j++ {
			// (*a.val)[index]
			value, _ := a.Getij(i, j)
			sum += value
		}
	}

	return
}

// trivial sum row/col
func Sum(a *Matrix, cursor int, vector uint8) (sum float64) {
	return Psum(a, cursor, vector, 0, 0, 0, 0)
}

// patrial sum row/col
func Psum(a *Matrix, cursor int, vector uint8, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	var n int

	if vector == ROWVECTOR {
		n = a.Cols() - hi_j
		for j := 0 + lo_j; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			sum += value
		}
	} else {
		n = a.Rows() - hi_i
		for i := 0 + lo_i; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			sum += value
		}
	}

	return
}

/*
SUMIF
Adds the cells specified by a given criteria. This function is used to browse a range when you search for a certain value.
The search supports regular expressions. You can enter "all.*", for example to find the first location of "all" followed by any characters. If you want to search for a text that is also a regular expression, you must precede every character with a \ character. You can switch the automatic evaluation of regular expression on and off in Tools - Options - LibreOffice Calc - Calculate.

Syntax
SUMIF(Range; Criteria; SumRange)

Range is the range to which the criteria are to be applied.
Criteria is the cell in which the search criterion is shown, or the search criterion itself. If the criteria is written into the formula, it has to be surrounded by double quotes.
SumRange is the range from which values are summed. If this parameter has not been indicated, the values found in the Range are summed.

SUMIF supports the reference concatenation operator (~) only in the Criteria parameter, and only if the optional SumRange parameter is not given.

Example
To sum up only negative numbers: =SUMIF(A1:A10;"<0")
=SUMIF(A1:A10;">0";B1:10) - sums values from the range B1:B10 only if the corresponding values in the range A1:A10 are >0.
See COUNTIF() for some more syntax examples that can be used with SUMIF().
*/

/*
SUMSQ
If you want to calculate the sum of the squares of numbers (totaling up of the squares of the arguments), enter these into the text fields.

Syntax
SUMSQ(Number1; Number2; ...; Number30)

Number1 to 30 are up to 30 arguments the sum of whose squares is to be calculated.

Example
If you enter the numbers 2; 3 and 4 in the Number 1; 2 and 3 text boxes, 29 is returned as the result.
*/

// sumsq of all table cells
func SumsqM(a *Matrix) (sum float64) {
	return PsumsqM(a, 0, 0, 0, 0)
}

// partial sumsq table
func PsumsqM(a *Matrix, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	for i := 0 + lo_i; i < a.Rows()-hi_i; i++ {
		for j := 0 + lo_j; j < a.Cols()-hi_j; j++ {
			value, _ := a.Getij(i, j)
			sum += value * value
		}
	}

	return
}

// trivial sumsq row/col
func Sumsq(a *Matrix, cursor int, vector uint8) (sum float64) {
	return Psumsq(a, cursor, vector, 0, 0, 0, 0)
}

// patrial sumsq row/col
func Psumsq(a *Matrix, cursor int, vector uint8, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	var n int

	if vector == ROWVECTOR {
		n = a.Cols() - hi_j
		for j := 0 + lo_j; j < n; j++ {
			value, _ := a.Getij(cursor, j)
			sum += value * value
		}
	} else {
		n = a.Rows() - hi_i
		for i := 0 + lo_i; i < n; i++ {
			value, _ := a.Getij(i, cursor)
			sum += value * value
		}
	}

	return
}

/*
SUMPRODUCT
Multiplies corresponding elements in the given arrays, and returns the sum of those products.

Syntax
SUMPRODUCT(Array1; Array2...Array30)

Array1, Array2...Array30 represent arrays whose corresponding elements are to be multiplied.
At least one array must be part of the argument list. If only one array is given, all array elements are summed.

Example

  A  B  C  D
1 2  3  4  5
2 6  7  8  9
3 10 11 12 13

=SUMPRODUCT(A1:B3;C1:D3) returns 397.
Calculation: A1*C1 + B1*D1 + A2*C2 + B2*D2 + A3*C3 + B3*D3
You can use SUMPRODUCT to calculate the scalar product of two vectors.

SUMPRODUCT returns a single number, it is not necessary to enter the function as an array function.
*/

func SumProduct(tables ...*Matrix) (sum float64) {
	// len(tables) > 1
	// If only one array is given, all array elements are summed.
	// all tables.Rows() are equals
	// all tables.Cols() are equals
	//for _, table := range tables {
	//	Mprintf("\ntable", "%.2f ", table)
	//}

	for i := 0; i < tables[0].Rows(); i++ {
		for j := 0; j < tables[0].Cols(); j++ {

			var prod float64 = 1

			for t := 0; t < len(tables); t++ {
				value, _ := tables[t].Getij(i, j)
				prod *= value
			}

			sum += prod
		}
	}

	return
}

/*
SUMX2MY2
Returns the sum of the difference of squares of corresponding values in two arrays.

Syntax
SUMX2MY2(ArrayX; ArrayY)

ArrayX represents the first array whose elements are to be squared and added.
ArrayY represents the second array whose elements are to be squared and subtracted.
*/
// разность квадратов
// sumx2my2 of all table cells
func Sumx2my2M(a, b *Matrix) (sum float64) {
	return Psumx2my2M(a, b, 0, 0, 0, 0)
}

// partial sumx2my2 table
func Psumx2my2M(a, b *Matrix, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	// a.Cols() == b.Cols()
	// a.Rows() == b.Rows()
	for i := 0 + lo_i; i < a.Rows()-hi_i; i++ {
		for j := 0 + lo_j; j < a.Cols()-hi_j; j++ {
			value_a, _ := a.Getij(i, j)
			value_b, _ := b.Getij(i, j)
			sum += value_a*value_a - value_b*value_b
		}
	}

	return
}

// trivial sumx2my2 row/col
func Sumx2my2(a, b *Matrix, cursor_a, cursor_b int, vector uint8) (sum float64) {
	return Psumx2my2(a, b, cursor_a, cursor_b, vector, 0, 0, 0, 0)
}

// patrial sumx2my2 row/col
func Psumx2my2(a, b *Matrix, cursor_a, cursor_b int, vector uint8, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	// a.Cols() == b.Cols()
	// a.Rows() == b.Rows()

	var n int

	if vector == ROWVECTOR {
		n = a.Cols() - hi_j
		for j := 0 + lo_j; j < n; j++ {
			value_a, _ := a.Getij(cursor_a, j)
			value_b, _ := a.Getij(cursor_b, j)
			sum += value_a*value_a - value_b*value_b
		}
	} else {
		n = a.Rows() - hi_i
		for i := 0 + lo_i; i < n; i++ {
			value_a, _ := a.Getij(i, cursor_a)
			value_b, _ := a.Getij(i, cursor_b)
			sum += value_a*value_a - value_b*value_b
		}
	}

	return
}

/*
SUMX2PY2
Returns the sum of the sum of squares of corresponding values in two arrays.

Syntax
SUMX2PY2(ArrayX; ArrayY)

ArrayX represents the first array whose elements are to be squared and added.
ArrayY represents the second array, whose elements are to be squared and added.
*/
// сумма квадратов
// sumx2py2 of all table cells
func Sumx2py2M(a, b *Matrix) (sum float64) {
	return Psumx2py2M(a, b, 0, 0, 0, 0)
}

// partial sumx2py2 table
func Psumx2py2M(a, b *Matrix, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	// a.Cols() == b.Cols()
	// a.Rows() == b.Rows()
	for i := 0 + lo_i; i < a.Rows()-hi_i; i++ {
		for j := 0 + lo_j; j < a.Cols()-hi_j; j++ {
			value_a, _ := a.Getij(i, j)
			value_b, _ := b.Getij(i, j)
			sum += value_a*value_a + value_b*value_b
		}
	}

	return
}

// trivial sumx2py2 row/col
func Sumx2py2(a, b *Matrix, cursor_a, cursor_b int, vector uint8) (sum float64) {
	return Psumx2py2(a, b, cursor_a, cursor_b, vector, 0, 0, 0, 0)
}

// patrial sumx2py2 row/col
func Psumx2py2(a, b *Matrix, cursor_a, cursor_b int, vector uint8, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	// a.Cols() == b.Cols()
	// a.Rows() == b.Rows()

	var n int

	if vector == ROWVECTOR {
		n = a.Cols() - hi_j
		for j := 0 + lo_j; j < n; j++ {
			value_a, _ := a.Getij(cursor_a, j)
			value_b, _ := a.Getij(cursor_b, j)
			sum += value_a*value_a + value_b*value_b
		}
	} else {
		n = a.Rows() - hi_i
		for i := 0 + lo_i; i < n; i++ {
			value_a, _ := a.Getij(i, cursor_a)
			value_b, _ := a.Getij(i, cursor_b)
			sum += value_a*value_a + value_b*value_b
		}
	}

	return
}

/*
SUMXMY2
Adds the squares of the variance between corresponding values in two arrays.

Syntax
SUMXMY2(ArrayX; ArrayY)

ArrayX represents the first array whose elements are to be subtracted and squared.
ArrayY represents the second array, whose elements are to be subtracted and squared.
*/
// квадрат разности

// sumxmy2 of all table cells
func Sumxmy2M(a, b *Matrix) (sum float64) {
	return Psumxmy2M(a, b, 0, 0, 0, 0)
}

// partial sumxmy2 table
func Psumxmy2M(a, b *Matrix, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	// a.Cols() == b.Cols()
	// a.Rows() == b.Rows()
	for i := 0 + lo_i; i < a.Rows()-hi_i; i++ {
		for j := 0 + lo_j; j < a.Cols()-hi_j; j++ {
			value_a, _ := a.Getij(i, j)
			value_b, _ := b.Getij(i, j)
			sum += (value_a - value_b) * (value_a - value_b)
		}
	}

	return
}

// trivial sumxmy2 row/col
func Sumxmy2(a, b *Matrix, cursor_a, cursor_b int, vector uint8) (sum float64) {
	return Psumxmy2(a, b, cursor_a, cursor_b, vector, 0, 0, 0, 0)
}

// patrial sumxmy2 row/col
func Psumxmy2(a, b *Matrix, cursor_a, cursor_b int, vector uint8, lo_i, lo_j, hi_i, hi_j int) (sum float64) {
	// a.Cols() == b.Cols()
	// a.Rows() == b.Rows()

	var n int

	if vector == ROWVECTOR {
		n = a.Cols() - hi_j
		for j := 0 + lo_j; j < n; j++ {
			value_a, _ := a.Getij(cursor_a, j)
			value_b, _ := a.Getij(cursor_b, j)
			sum += (value_a - value_b) * (value_a - value_b)
		}
	} else {
		n = a.Rows() - hi_i
		for i := 0 + lo_i; i < n; i++ {
			value_a, _ := a.Getij(i, cursor_a)
			value_b, _ := a.Getij(i, cursor_b)
			sum += (value_a - value_b) * (value_a - value_b)
		}
	}

	return
}
