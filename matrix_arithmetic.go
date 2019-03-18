//
// Pure Go implementation of the Matrix arithmetic and GSL libs
// draft version
// unstable API
//
// Linear algebra(base)
//
// If you need just good matrix lib you may use this one 
// https://godoc.org/github.com/wiless/go.matrix
// This solution consists of little more other things than 
// just matrices and determinants
//
package matrix_arithmetic

import(
	// common purpose
	"fmt"
	"bufio"
	_ "os"
	"io"
	"strings"
	"strconv" // aoti and so on convertions
	"errors" // errors.New()
	"math"
	_ "log"
	_ "time"
//	"sync"
)

const(

	zero = 1.0e-7 // zero tolerance


	LINELEN = 4096

	/* VECTOR TYPES */
	ROWVECTOR = 0
	COLVECTOR = 1

	/* ERRORS */
	//TODO: fmt.Errorf("Value not found for type %v", argType
	RMISMATCH = "row mismatch"
	CMISMATCH = "column mismatch"
	NOTSQUARE = "not a square Matrix"
	ALLOCFAIL = "allocation failure"
	PARSEFAIL = "parse file failure"
	FILEREADFAIL = "file read failure"
	ROWPARSEFAIL = "row parse failure"
	COLPARSEFAIL = "column parse failure"
	RCMISMATCH = "row-column mismatch"
	INDEXOUTOFRANGE = "index out of range"
	LENMISMATCH = "length mismatch"
	NULLARG = "NULL argument"
)

var(

)

/*
 * Matrix interface
 * export 
 */
type IMatrix interface {
	String() string
	Rows() int
	Cols() int
	Getij(int, int) (float64, error)
	Setij(int, int, float64) error
	//TODO: m_ops
	//<...>

	// interface for GNU GSL Statistics library (v1.15, GPLv3) implemented in Go
	Get(int) float64
	Len() int
}

/*
 * Matrix structure 
 */
type Matrix struct {
	rows int
	cols int
	val *[]float64
}

/*
 * Matrix methods
 */
func (m *Matrix) Rows() int {
	return m.rows
}

func (m *Matrix) Cols() int {
	return m.cols
}

func (m *Matrix) Getij(i int, j int) (float64, error) {
	return Getij(m, i, j)
}

func (m *Matrix) Setij(i int, j int, value float64) error {
	return Setij(m, i, j, value)
}


/*TODO: ??? interface for GNU GSL Statistics library (v1.15, GPLv3) implemented in Go */
func (m *Matrix) Get(index int) float64 {
	return float64(0)
}
func (m *Matrix) Len() int {
	return 0
}



/* allocates a new Matrix, elements not initialized */
// export
func New(nrows, ncols int) (*Matrix, error) {
	if nrows <= 0 || ncols <= 0 {
		return nil,errors.New(ALLOCFAIL)
	}

	// create Matrix struct
	var m *Matrix = new(Matrix)

	// init Matrix
	m.rows = nrows
	m.cols = ncols
	m.val = new([]float64)

	// create nrows x ncols region
	*(m.val) = make([]float64, m.rows * m.cols )

	// success
	return m, nil
}


/* allocates a new Matrix and initializes it from a file */
// export
func Fnew(source io.Reader) ([]*Matrix, error) {
	//TODO: this oldschool matrices gethering must be rewieved in the future

	//TODO: vebose parseing

	reader := bufio.NewReader(source)

	var lineno uint = 0
	var rows int = 0
	var cols int = 0
	var index int
	var matrices []*Matrix
	var matrix = -1 // just for starting from zero at the first time
	stop := false

	for ;!stop; {
		switch line, err := reader.ReadString('\n'); (err) {
		case io.EOF: // end-of-file
			stop = true
		case nil: // success
			lineno++

			splitted := strings.Split(line[:len(line)-1], ",")

			// skip emty line
			if len(splitted[0]) <= 0 {
				fmt.Println("SKIP EMPTY LINE")
				goto NEXT_LINE
			}

			// skip comment
			if splitted[0][0] == '#' {
				fmt.Println("SKIP COMMENT:", splitted)
				goto NEXT_LINE
			}

			// get number of rows
			if splitted[0] == "rows" {
				i, err := strconv.Atoi(splitted[1])
				if err != nil {
					//handle convetion error
					fmt.Println(err)
					return nil,err
				}
				rows = i
				fmt.Println("RAWS", rows)
				goto NEXT_LINE
			}

			// get number of columns
			if splitted[0] == "cols" {
				i, err := strconv.Atoi(splitted[1])
				if err != nil {
					//handle convetion error
					fmt.Println(err, lineno)
					//fmt.Errorf("Value not found for type %v", argType
					return nil,err
				}
				cols = i
				fmt.Println("COLS", cols)

				// append to creating Matrixes list
				if rows > 0 && cols > 0 {
					// reset index at the beginnig
					index = 0
					// allocate new Matrix
					fmt.Printf("ALLOCATE NEW MATRIX %d x %d\n", rows, cols)
					a, err := New(rows, cols)
					if err != nil {
						fmt.Println(err) // error checking done in m_new
						return nil,err
					}
					// add new mantrix
					matrices = append(matrices, a)
					// setup poiter
					matrix++
					goto NEXT_LINE
				} else {
					// handle length error
					fmt.Println("INDEXOUTOFRANGE", lineno)
					return nil,errors.New(INDEXOUTOFRANGE)
				}
			}

			// fill row's cells of the Matrix
			if len(splitted) == cols {
				for j := 0; j < cols; j++ {
					// covert to float64
					value, err := strconv.ParseFloat(splitted[j], 64)
					if err != nil {
						//handle convetion error
						fmt.Println("CONVFAIL", err, lineno)
						return nil,err
					}
					// insert into Matrix
					(*matrices[matrix].val)[index] = value
					index++
				}
				//fmt.Println()
			} else {
				fmt.Println("PARSEFAIL line:", lineno)
				stop = true
				return nil,errors.New(PARSEFAIL)
			}
NEXT_LINE:
			break
		default: // error 
			fmt.Println("FILEREADFAIL", err, lineno)
			stop = true
			return nil,errors.New(FILEREADFAIL)
		}
	}

	return matrices,nil
}


/* Matrix comparison a == b */
// TODO maybe add random control method too? 
// export
func Compare(a *Matrix, b *Matrix) (bool,error) {
	var result bool = false

	if a == nil || b == nil {
		return result, errors.New(NULLARG)
	}
	if a.rows != b.rows {
		return result, errors.New(RMISMATCH)
	}
	if a.cols != b.cols {
		return result, errors.New(CMISMATCH)
	}

	// matrices are equal apriori
	result = true
	stop := false
	for i := 0; i < a.rows * a.cols && !stop; i++ {
		if (*a.val)[i] != (*b.val)[i] {
			stop = true
			result = false
		}
	}

	// success
	return result, nil
}


/* Matrix addition sum = a + b */
// export
func Add(sum *Matrix, a *Matrix, b *Matrix) error {
	if sum == nil || a == nil || b == nil {
		return errors.New(NULLARG)
	}
	if a.rows != b.rows {
		return errors.New(RMISMATCH)
	}
	if sum.rows != b.rows {
		return errors.New(RMISMATCH)
	}
	if a.cols != b.cols {
		return errors.New(CMISMATCH)
	}
	if b.cols != sum.cols {
		return errors.New(CMISMATCH)
	}

	index := 0
	for i := 0; i < sum.rows; i++ {
		for j := 0; j < sum.cols; j++ {
			(*sum.val)[index] = (*a.val)[index] + (*b.val)[index]
			index++
		}
	}

	// success
	return nil
}


/* Matrix subtraction diff = a - b */
// export
func Sub(diff *Matrix, a *Matrix, b *Matrix) error {
	if diff == nil || a == nil || b == nil {
		return errors.New(NULLARG)
	}
	if a.rows != b.rows {
		return errors.New(RMISMATCH)
	}
	if diff.rows != b.rows {
		return errors.New(RMISMATCH)
	}
	if a.cols != b.cols {
		return errors.New(CMISMATCH)
	}
	if diff.cols != b.cols {
		return errors.New(CMISMATCH)
	}

	index := 0
	for i := 0; i < diff.rows; i++ {
		for j := 0; j < diff.cols; j++ {
			(*diff.val)[index] = (*a.val)[index] - (*b.val)[index]
			index++
		}
	}

	// success
	return nil
}


/* Matrix assignment a = b */
// export
func Assign(a *Matrix, b *Matrix) error {
	if a == nil || b == nil {
		return errors.New(NULLARG)
	}
	if a.rows != b.rows {
		return errors.New(RMISMATCH)
	}
	if a.cols != b.cols {
		return errors.New(CMISMATCH)
	}

	index := 0
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			(*a.val)[index] = (*b.val)[index]
			index++
		}
	}

	// success
	return nil
}


/* Matrix multiplication prod = a * b */
// export
func Mult(prod *Matrix, a *Matrix, b *Matrix) error {
	if prod == nil || a == nil || b == nil {
		return errors.New(NULLARG)
	}
	if prod.rows != a.rows {
		return errors.New(RMISMATCH)
	}
	if prod.cols != b.cols {
		return errors.New(CMISMATCH)
	}
	if a.cols != b.rows {
		return errors.New(RMISMATCH)
	}

	for i := 0; i < a.rows; i++ {
		for j := 0; j < b.cols; j++ {
			(*prod.val)[mdx(prod,i,j)] = float64(0)
			for k := 0; k < a.cols; k++ {
				(*prod.val)[mdx(prod,i,j)] += (*a.val)[mdx(a,i,k)] * (*b.val)[mdx(b,k,j)]
			}
		}
	}

	// success
	return nil

}


/* Matrix transposition trans = a transpose */
// export
func Transpose(trans *Matrix, a *Matrix) error {
	if trans == nil || a == nil {
		return errors.New(NULLARG)
	}
	if trans.rows != a.cols {
		return errors.New(RMISMATCH)
	}
	if trans.cols != a.rows {
		return errors.New(RMISMATCH)
	}

	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			(*trans.val)[mdx(trans,j,i)] = (*a.val)[mdx(a,i,j)]
		}
	}

	// success
	return nil
}


func Tprintf(comment string, nformat string, a *Matrix, b *Matrix) error {
	if a == nil || b == nil {
		return errors.New(NULLARG)
	}

	fmt.Printf("%s\n", comment)

	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			fmt.Printf(nformat, (*a.val)[mdx(a,i,j)])
		}
		fmt.Printf(" ** ");
		for j := 0; j < a.cols; j++ {
			fmt.Printf(nformat, (*b.val)[mdx(a,i,j)])
		}
		fmt.Println()
	}
	fmt.Println()

	// success
	return nil
}


/* assigns the values of the identity Matrix to iden */
// export
func AssignIdentity(iden *Matrix) error { // or maybe (iden *Matrix, result error) ?
	if iden == nil {
		return errors.New(NULLARG)
	}
	if iden.rows != iden.cols {
		return errors.New(NOTSQUARE)
	}

	index := 0;
	for i := 0; i < iden.rows; i++ {
		for j := 0; j < iden.cols; j++ {
			(*iden.val)[index] = 0.0;
			if (i == j) {
				(*iden.val)[index] = 1.0
			}
			index++;
		}
	}

	// success
	return nil;
}


/* returns the value of a specified Matrix element */
// export
func Getij(a *Matrix, i int, j int) (float64,error) {
	if a == nil {
		return float64(0), errors.New(NULLARG)
	}
	if i > a.rows-1 {
		return float64(0), errors.New(INDEXOUTOFRANGE)
	}
	if j > a.cols-1 {
		return float64(0), errors.New(INDEXOUTOFRANGE)
	}

	// success
	return (*a.val)[mdx(a, i, j)], nil
}


/* sets a specified Matrix element to specified value */
// export
func Setij(a *Matrix, i int, j int, value float64) error {
	if a == nil {
		return errors.New(NULLARG)
	}
	if i > a.rows-1 {
		return errors.New(INDEXOUTOFRANGE)
	}
	if j > a.cols-1 {
		return errors.New(INDEXOUTOFRANGE)
	}

	(*a.val)[mdx(a, i, j)] = value

	// success
	return nil
}


/* returns the number of rows in a Matrix */
// get_rows_count int?
// export
func GetRows(a *Matrix) (int,error) {
	if a == nil {
		return 0, errors.New(NULLARG)
	}

	// success
	return a.rows, nil
}


/* returns the number of columns in a Matrix */
// get_cols_count int ?
// export
func GetCols(a *Matrix) (int,error) {
	if a == nil {
		return 0, errors.New(NULLARG)
	}
	return a.cols, nil
}


/* prints a Matrix to stdout with a label using a specified numeric format */
// export
func Mprintf(label string, format string, a *Matrix) error {
	if a == nil {
		return errors.New(NULLARG)
	}
	if format == "" {
		format = "[%f]"
	}

	fmt.Printf("%s\n", label);
	fmt.Printf("rows = %d, cols = %d\n", a.rows, a.cols)
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			fmt.Printf(format, (*a.val)[mdx(a,i,j)])
		}
		fmt.Println();
	}

	// success
	return nil
}


/* write a Matrix a into a fp file in a CSV format */
// export
func Fputcsv(fp io.Writer, a *Matrix) error {
	var sep string
	var comma string = ","
	var nocomma string = ""

	if a == nil {
		return errors.New(NULLARG)
	}

	fmt.Fprintf(fp, "rows,%d\n", a.rows);
	fmt.Fprintf(fp, "cols,%d\n", a.cols);
	for i := 0; i < a.rows; i++ {
		sep = nocomma;
		for j := 0; j < a.cols; j++ {
			fmt.Fprintf(fp, "%s%f", sep, (*a.val)[mdx(a,i,j)]);
			sep = comma;
		}
		fmt.Fprintf(fp, "\n");
	}

	// success
	return nil
}


/* assigns the element values of a Matrix from a 2D array */
// export
/*TODO: ???
func AssignArr2D(a *Matrix, nrows int, ncols int, arr *[]float64) error {
	if a == nil {
		return errors.New(NULLARG)
	}
	if nrows != a.rows {
		return errors.New(RMISMATCH)
	}
	if ncols != a.cols {
		return errors.New(CMISMATCH)
	}

	length := nrows * ncols

	return nil //AssignArr1D(a, length, arr)
}*/
func AssignArray2D(a *Matrix, nrows int, ncols int, arr *[][]float64) error {
	if a == nil {
		return errors.New(NULLARG)
	}
	if nrows != a.rows {
		return errors.New(RMISMATCH)
	}
	if ncols != a.cols {
		return errors.New(CMISMATCH)
	}

	index := 0
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			(*a.val)[index] = (*arr)[i][j]
			index++
		}
	}

	// success
	return nil
}

/* assigns the element values of a Matrix from a 1D array */
// export
func AssignArray1D(a *Matrix, alen int, arr *[]float64) error {
	if a == nil {
		return errors.New(NULLARG)
	}
	if alen != (a.rows * a.cols) {
		return errors.New(LENMISMATCH)
	}

	index := 0
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			(*a.val)[index] = (*arr)[index]
			index++
		}
	}

	// success
	return nil
}


/* multiply each element of Matrix by a constant c 
	 prod = c * a */
// export
func MultConst(cprod *Matrix, c float64, a *Matrix) error {
	if cprod == nil || a == nil {
		return errors.New(NULLARG)
	}
	if cprod.rows != a.rows {
		return errors.New(RMISMATCH)
	}
	if cprod.cols != a.cols {
		return errors.New(CMISMATCH)
	}

	index := 0
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			(*cprod.val)[index] = c * (*a.val)[index]
			index++
		}
	}

	// success
	return nil
}


/* returns the absolute value of the Matrix element
   with the largest absolute value */
// export
func MaxAbsElement(a *Matrix) (float64, error) {
	var maxelement float64 = float64(0)
	if a == nil {
		return float64(0),errors.New(NULLARG)
	}

	index := 0
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			if math.Abs((*a.val)[index]) > maxelement {
				maxelement = math.Abs((*a.val)[index])
			}
			index++
		}
	}

	// success
	return maxelement, nil
}


/* returns the trace of the Matrix */
// export
func Trace(a *Matrix) (float64, error) {
	var tr = float64(0)
	if a == nil {
		return float64(0), errors.New(NULLARG)
	}

	for i := 0; i < a.rows; i++ {
		j := i
		tr += (*a.val)[mdx(a,i,j)]
	}

	// sucecss
	return tr, nil
}


/* returns the Euclidean norm of the Matrix, i.e.,
   the square root of the sum of the squares of
   each element */
// export
func E_Norm(a *Matrix) (float64, error) {
	var norm = float64(0)
	if a == nil {
		return float64(0), errors.New(NULLARG)
	}

	index := 0
	for i := 0; i < a.rows; i++ {
		for j := 0; j < a.cols; j++ {
			norm += (*a.val)[index] * (*a.val)[index]
			index++
		}
	}

	// sucecss
	return math.Sqrt(norm), nil
}


/* Matrix inversion inv = a inverse 
 * returns the value of determinat */
// export
func Inverse(v *Matrix, a *Matrix, epsilon float64) (float64,error) {
	/* calculates the inverse of Matrix a using Gauss-Jordan
	   elimination with partial pivot maximization */
	var col int
	var row int
	var pivot float64
	var swap int
	var det float64
	var sign int

	if v == nil || a == nil {
		return float64(0), errors.New(NULLARG)
	}
	if v.rows != a.rows {
		return float64(0), errors.New(RMISMATCH)
	}
	if v.cols != a.cols {
		return float64(0), errors.New(CMISMATCH)
	}
	if v.rows != v.cols {
		return float64(0), errors.New(NOTSQUARE)
	}

	// initialize
	e_norm,_ := E_Norm(a);
	det = float64(1);
	sign = 1;

	/* allocate a "scratch" Matrix to invert */
	t, err := New(a.rows, a.cols)
	if err != nil {
		return float64(0), err
	}

	Assign(t, a);

	/* set target Matrix to the identity Matrix */
	AssignIdentity(v);

	for row = 0; row < t.rows; row++ {
		/* find largest element below diagonal in column */
		swap = maxelementrow(t, row)
		/* swap rows to put largest element on pivot */
		if swap != row {
			if sign > 0 {
				sign = -1
			} else {
				sign = 1
			}
			swaprows2(t, v, row, swap)
		}
		/* divide each element on pivot row by pivot
		 * element putting a 1 in the pivot element */
		pivot = (*t.val)[mdx(t, row, row)]
		det *= pivot
		if (math.Abs(det) / e_norm) < epsilon {
			return float64(0), nil;  /* potentially singular Matrix */
		}

		multrow(t, row, float64(1) / pivot)
		multrow(v, row, float64(1) / pivot)

		/* subtract a multiple of the pivot row from each
		 * row to put 0's in all elements of the pivot
		 * column except the pivot row */
		col = row
		set_col_zero(t, v, col)
	}

	if sign < 0 {
		det = -det
	}

	// success
	return det, nil
}




/* returns the value of the determinant */
// export
func Det(a *Matrix, epsilon float64) (float64, error) {
	/* calculates the determinant of Matrix a using Gaussian
	 * elimination with partial pivot maximization */
	var pivot float64
	var e_norm float64
	var sign int
	var col int
	var row int
	var det float64

	if a == nil {
		return float64(0),errors.New(NULLARG)
	}
	if a.rows != a.cols {
		return float64(0),errors.New(NOTSQUARE)
	}

	// initialize
	e_norm, _ = E_Norm(a)
	det = float64(1)
	sign = 1

	/* allocate a "scratch" Matrix to work with */
	t, err := New(a.rows, a.cols)
	if err != nil {
		return float64(0), err
	}
	Assign(t, a)

	/* for each row */
	for row = 0; row < t.rows; row++ {
		/* find largest element below diagonal in column */
		swap := maxelementrow(t, row)
		/* swap rows to put largest element on pivot */
		if swap != row {
			sign = -sign
			swaprows(t, row, swap)
		}
		/* multiply running product of det by the pivot
		 * element */
		pivot = (*t.val)[mdx(t, row, row)]
		det *= pivot
		if (math.Abs(det) / e_norm) < epsilon {
			return float64(0), nil  /* potentially singular Matrix */
		}
		/* subtract a multiple of the pivot row from each
		 * row (below the diagonal) to put 0's in all
		 * elements of the pivot column below the diagonal */
		col = row
		set_low_zero(t, col)
	}

	if sign < 0 {
		det = -det
	}
	// success
	return det, nil
}



/**
 ** Internals
 */

/* function to compute index of 1D array corresponding to
 * i,j indices of 2D array 
 * mdx is defined as a function for development
 * and testing; it is commented out and replaced
 * with a macro (a.cols * i + j) for production */
func mdx(a *Matrix, i int, j int) int {
	return i * a.cols + j
}


/* subtracts multiples of pivot row from each row below the 
 * pivot to make each column element below the diagonal in
 * column "col" equal to zero */
func set_low_zero(t *Matrix, col int) {
	pivot := (*t.val)[mdx(t, col, col)]
	for i := col + 1; i < t.rows; i++ {
		factor := (*t.val)[mdx(t, i, col)] / pivot
		for j := col; j < t.cols; j++ {
			(*t.val)[mdx(t, i, j)] -= factor * (*t.val)[mdx(t, col, j)]
		}
	}
}


/* swaps rows "row" and "swap" in Matrix t */
func swaprows(t *Matrix, row int, swap int) {
	for j := 0; j < t.cols; j++ {
		temp := (*t.val)[mdx(t, row, j)];
		(*t.val)[mdx(t, row, j)] = (*t.val)[mdx(t, swap, j)];
		(*t.val)[mdx(t, swap, j)] = temp;
	}
}


/* function to find the row (on or below the diagonal) in
 * which the maximum element of a column occurs */
func maxelementrow(t *Matrix, col int) int {
	trial := float64(0)
	maximum := float64(0)
	row := col
	ilargest := row

	for i := row; i < t.rows; i++ {
		trial = math.Abs((*t.val)[mdx(t, i, col)])
		if trial > maximum {
			maximum = trial
			ilargest = i
		}
	}

	return ilargest
}


/* function to subtract multiples of the pivot row from all 
 * other rows to force each element in the pivot column to
 * * zero (except the actual pivot element) */
func set_col_zero(t *Matrix, v *Matrix, col int) {
	pivot := (*t.val)[mdx(t, col, col)]
	for i := 0; i < t.rows; i++ {
		if i == col {
			continue
		}
		factor := (*t.val)[mdx(t, i, col)] / pivot
		for j := 0; j < t.cols; j++ {
			(*t.val)[mdx(t, i, j)] -= factor * (*t.val)[mdx(t, col, j)];
			(*v.val)[mdx(v, i, j)] -= factor * (*v.val)[mdx(v, col, j)];
		}
	}
}


/* swap row numbers row and swap in matrices t and v */
func swaprows2(t *Matrix, v *Matrix, row int, swap int) {
	for j := 0; j < t.cols; j++ {
		temp := (*t.val)[mdx(t, row, j)]
		(*t.val)[mdx(t, row, j)] = (*t.val)[mdx(t, swap, j)]
		(*t.val)[mdx(t, swap, j)] = temp
		temp = (*v.val)[mdx(v, row, j)]
		(*v.val)[mdx(v, row, j)] = (*v.val)[mdx(v, swap, j)]
		(*v.val)[mdx(v, swap, j)] = temp
	}
}


/* multiplies the specified row in Matrix a by a constant * "c" */
func multrow(a *Matrix, row int, c float64) error {
	if a == nil {
		return errors.New(NULLARG)
	}
	index := row * a.cols;
	for j := 0; j < a.cols; j++ {
		(*a.val)[index] *= c;
		index++;
	}
	return nil;
}

