//
// Pure Go implementation of the Matrix arithmetic and GSL libs
// draft version
// unstable API
//
// Methods and algorithm drafts and sources
//
// Start point of algorithms and methods implementation.
// Examples of specific algorithms for test are over here too.
//
package matrix_arithmetic


import(
	// common purpose
	_"fmt"
	_ "bufio"
	_"os"
	_ "io"
	_ "strings"
	_ "strconv" // aoti and so on convertions
	_ "errors" // errors.New()
	_"math"
	_ "log"
	_"time"
//	"sync"
)
/*
// n-мерная норма
// условие остановки
func norma(cur, next []float64) float64 {
	var sumSquares float64

	for i := 0; i < len(cur); i++ {
		sumSquares += (next[i] - cur[i])*(next[i] - cur[i])
	}

	//fmt.Println("    norma:", sumSquares)
	//fmt.Println("    norma:", math.Sqrt(sumSquares))
	return math.Sqrt(sumSquares)
}


// Метод половинного деления для нахождения минимума в градиентном спуске
// dichotomia — разделение на две части
func dichotomia(g func(args []float64, part float64) float64, args []float64, a0, b0, epsilon float64) float64 {
	// Номер шага
	var k int
	// Отклонени от середины отрезка влево, вправо
	var left float64
	var right float64
	// Величина на которую мы отклонимся от середины отрезка
	var deviation float64 = 0.5 * epsilon
	// Точка минимума
	var x_min float64
	// Отрезок локализации минимума
	var ak float64 = a0
	var bk float64 = b0

	//	Шаг 3. Вычислить grad f(x[k]).
	//	Шаг 4. Проверить выполнение критерия окончания  ||grad f(x[k])|| < epsilon1 :
	//		а) если критерий выполнен, то x[0] = x[k], останов;
	//		б) если критерий не выполнен, то перейти к шагу 5.
	// Пока длина отрезка больше заданной точности
	for k = 1; (bk - ak) >= epsilon; k++ {
		// Берем середину (ну почти середину - +\- некоторое дельта 
		// в частности у нас deviation = 0.5 * epsilon)
		left = (ak + bk - deviation) / 2.0
		right = (ak + bk + deviation) / 2.0

		// Шаг 3. Вычислить grad f(x[k])
		// Проверяем в какую часть попадает точка минимума для этого 
		// вычисляем функцию g для точек слева и справа от разбиения 
		// сравния их определяем полжение точки, слева или справа
		// и выбираем соответствующую найденную точку
		if g(args, left) <= g(args, right) {
		// Теперь правая граница отрезка локализации равна right
			bk = right;
		} else {
		// Теперь левая граница отрезка локализации равна left
			ak = left;
		}
	}

	//делим получившийся отрезок пополам и получаем точку минимума 
	x_min = (ak + bk) / 2.0;

	return x_min;
}


func grad(function func (args []float64) float64, args []float64, i int, delta float64) float64 {
	// NOTE: delta ought to be small enough but you should remember 
	//       that too small value will drive to reducing accuracy
	//
	// df/dxi = f(x1,x2,...,xi+/\xi,...xn) - f(x1,x2,...xi-/\xi,...xn) / (2 * /\xi) 
	//
		left := make([]float64, len(args))
		copy(left, args)
		right := make([]float64, len(args))
		copy(right, args)
		left[i] += delta
		right[i] -= delta
		return ( function(left) - function(right) ) / (2.0 * delta)
}

// Метод наискорейшего спуска
// steepest descent method
func SteepestDescent(function func(args []float64) float64, args []float64, epsilon float64) (float64,[]float64) {
	const MAXITERATIONS = 1000
	// приближение 
	var k int
	// вечина шага приближения
	var lambda float64

	// П.1. Задают начальное приближение и точность расчёта vec x[0], epsilon
	//  Шаг 1. Задать x[0], epsilon1 > 0, epsilon2 > 0, предельное число итераций М. 
	//   Найти градиент функции в произвольной точке. Определить частные производные функции f(x):
	//   grad f(x) = [ df(x)/x1,...,df(x)/dxn ]T(транспонированный вектор)
	// Начальное приближение u[0]
	u_cur := make([]float64, len(args))
	// Новое прилижение u[k]
	u_next := make([]float64, len(u_cur))

	// задать начальное приближение
	copy(u_cur, args)

	fmt.Println("## Исходная точка:")
	for i,x := range u_cur {
		fmt.Printf(" %d:  x[%d]:%.2f, ", k, i, x)
	}
	fmt.Println()

	fmt.Println()
	fmt.Println("## Приближения:")

	// Шаг 2. Положить k = 0
	stop := false
	for k = 0; k < MAXITERATIONS && !stop ; k++ {
		// П.2. Рассчитывают vec x[k+1] = vec x[k] - lambda[k] * nabla F(vec x[k]),
		// где lambda[k]= argmin_lambda F(vec x[k] - lambda * nabla F(vec x[k]))
		//   Шаг 3. Вычислить grad f(x[k]).
		//
		//   Шаг 4. Проверить выполнение критерия окончания  ||grad f(x[k])|| < epsilon1 :
		//   а) если критерий выполнен, то x[0] = x[k], останов;
		//   б) если критерий не выполнен, то перейти к шагу 5.
		//
		//   Шаг 5. Проверить выполнение неравенства k ≥ M:
		//   а) если неравенство выполнено, то x[0] = x[k], останов;
		//   б) если нет, то перейти к шагу 6
		//
		//   Шаг 6. Вычислить величину шага lambda[k0](на начальном шаге, k = 0) из условия
		//   F(lambda[k]) = f(x[k] - lambda[k] * grad f(x[k])) -> min lambda[k]
		//
		// Выберем метод дихотомии как F(lambda[k])
		argmin := dichotomia
		g := func (args []float64, lambda float64) float64 {
			// Это функция g(x[k]) в методе наискорейшего (градиентного) спуска
			// применяемая для нахождения шага lambda[k] как минимума функции g(x[k])
			// на отрезке(здес a=-10000, b=100000, методом дихотомии)

			temp := make([]float64, len(args))
			// шаг дифференцирования
			delta := 0.05 //TODO: ajuset
			for i := 0; i < len(args); i++ {
				temp[i] = args[i] - lambda * grad(function, args, i, delta)
			}

			return function(temp)
		}

		// Находим lambda[k] как минимум функции g(x[k]) на отрезке -10000,100000
		lambda = argmin(g, u_cur, -10000, 100000, epsilon);
		// Вычисляем u[k] новое прилижение
		// Шаг 7. Вычислить x[k+1] = x[k] - lambda[k] * grad f(x[k])
		// шаг дифференцирования
		delta := 0.05 //TODO: ajuset
		fmt.Printf(" %d: ", k)
		for i := 0; i < len(u_cur); i++ {
			u_next[i] = u_cur[i] - lambda * grad(function, u_cur, i, delta)
			fmt.Printf(" x[%d]: %.2f, ", i, u_next[i] )
		}
		fmt.Println()

		// П.3. Проверяют условие остановки:
		// Если |vec x[k+1] - vec x[k]| > epsilon
		//  или |F(vec x[k+1]) - F(vec x[k])| > epsilon
		//  или |nabla F(vec x[k+1])| > epsilon (выбирают одно из условий), 
		// то k = k + 1 и переход к П.2.
		// Иначе vec x = vec x[k+1] и останов.
		//   Шаг 8. Проверить выполнение условий ||x[k+1] - x[k]|| < epsilon2, ||f(x[k+1]) - f(x[k])|| < epsilon2:
		//   а) если оба условия выполнены при текущем значении k и k = k - 1, то расчет окончен, x[0] = x[k+1], останов;
		//   б) если хотя бы одно из условий не выполнено, то положить k = k + 1 и перейти к шагу 3
		if k > 1 {
			// Проверяем условие остановки
			if norma(u_cur, u_next) < epsilon {
				// останов
				stop = true
			}
		}

		copy(u_cur, u_next)
	} // eof for k

	fmt.Println()
	fmt.Println("## Точка минимума epsilon:", epsilon)
	fmt.Println(u_cur)
	fmt.Printf("f(x): %f\n", function(u_cur))

	// получить минимум функции
	return function(u_cur),u_cur
}





// ======================================================================
// Пример для проверки 
// 

// Собственно здесь записывается наша функция
func f(x,y float64) float64 {
	//TODO: здеась можно задать исследуемею функцию или поиграться с этой
	//NOTE: если поменяете функцию не забудьте про частные производные
	return x*x + 2.0*x*y + 3.0*y*y - 2.0*x -3.0*y
}
func f2(args []float64) float64 {
	//fmt.Println("f2 args:", args, len(args))
	return args[0]*args[0] + 2.0*args[0]*args[1] + 3.0*args[1]*args[1] - 2.0*args[0] - 3.0*args[1]
}

//
// операция grad или nabla
//
// Это первая производная по dx
// частная производная по х
func f_dx(x, y float64) float64 {
    return 2.0*x + 2.0*y - 2.0
}
// Это первая производная по dy
// частная производная по y
func f_dy(x, y float64) float64 {
    return 2.0*x + 6.0*y - 3.0
}

// Это функция g(x[k]) в методе наискорейшего (градиентного) спуска
// применяемая для нахождения шага lambda[k] как минимума функции g(x[k])
// на отрезке(здес a=-10000, b=100000) (здесь методом дихотомии)
func gExample(x, y, lambda float64) float64 {
	//var model_f_dx float64 = f_dx(x, y)
	//var model_f_dy float64 = f_dy(x, y)
	// срединный разностный метод вычисления первой производной
	// numerical differentiation
	// NOTE: delta ought to be small enough but you should remember 
	//       that too small value will drive to reducing accuracy
	//var delta float64 = 0.5
	//var approx_f_dx float64 = ( f(x + delta / 2.0, y) - f(x - delta / 2.0, y) ) / delta
	//var approx_f_dy float64 = ( f(x, y + delta / 2.0) - f(x, y - delta / 2.0) ) / delta

	//var args []float64 = []float64{x, y}
	//fmt.Println(" args before:", args)
	//var grad_dx = grad(f2, args, 0, 0.05)
	//var grad_dy = grad(f2, args, 1, 0.05)

	//fmt.Printf("  model x: %.5f  approx1 x: %.5f  grad x: %.5f\n", model_f_dx, approx_f_dx, grad_dx)
	//fmt.Printf("  model y: %.5f  approx1 y: %.5f  grad y: %.5f\n", model_f_dy, approx_f_dy, grad_dy)
	//fmt.Println()


	// Проверка дифференцирования    
	return f(x - lambda * f_dx(x, y), y - lambda * f_dy(x, y))
	//return f(x - lambda * approx_f_dx, y - lambda * approx_f_dy)
	//return f(x - lambda * grad(f2, args, 0, 0.05), y - lambda * grad(f2, args, 1, 0.05))
}


// двумерная норма
// условие остановки
func normaExample(x, y float64) float64 {
//	fmt.Println("    normaE:", x*x + y*y)
//	fmt.Println("    normaE:", math.Sqrt(x*x + y*y))
	return math.Sqrt(x*x + y*y)
}


// Метод половинного деления для нахождения минимума в градиентном спуске
// dichotomia — разделение на две части
func dichotomiaExample(a0, b0, epsilon, x, y float64) float64 {
	// Номер шага
	var k int
	// Отклонени от середины отрезка влево, вправо
	var left float64
	var right float64
	// Величина на которую мы отклонимся от середины отрезка
	var deviation float64 = 0.5 * epsilon
	// Точка минимума
	var x_min float64
	// Отрезок локализации минимума
	var ak float64 = a0
	var bk float64 = b0

	//	Шаг 3. Вычислить grad f(x[k]).
	//	Шаг 4. Проверить выполнение критерия окончания  ||grad f(x[k])|| < epsilon1 :
	//		а) если критерий выполнен, то x[0] = x[k], останов;
	//		б) если критерий не выполнен, то перейти к шагу 5.
	// Пока длина отрезка больше заданной точности
	for k = 1; (bk - ak) >= epsilon; k++ {
		// Берем середину (ну почти середину - +\- некоторое дельта 
		// в частности у нас deviation = 0.5 * epsilon)
		left = (ak + bk - deviation) / 2.0
		right = (ak + bk + deviation) / 2.0

		// Шаг 3. Вычислить grad f(x[k])
		// Проверяем в какую часть попадает точка минимума для этого 
		// вычисляем функцию g для точек слева и справа от разбиения 
		// сравния их определяем полжение точки, слева или справа
		// и выбираем соответствующую найденную точку
		if gExample(x, y, left) <= gExample(x, y, right) {
		// Теперь правая граница отрезка локализации равна right
			bk = right;
		} else {
		// Теперь левая граница отрезка локализации равна left
			ak = left;
		}
	}

	//делим получившийся отрезок пополам и получаем точку минимума 
	x_min = (ak + bk) / 2.0;

	return x_min;
}

// Метод наискорейшего спуска(Пример)
// steepest descent method
func SteepestDescentExample(bx, by, epsilon float64) float64 {
	const MAXITERATIONS = 1000
	var k int
	var x_cur float64
	var y_cur float64
	var x_next float64
	var y_next float64
	var lambda float64

	// П.1. Задают начальное приближение и точность расчёта vec x[0], epsilon
	//  Шаг 1. Задать x[0], epsilon1 > 0, epsilon2 > 0, предельное число итераций М. 
	//   Найти градиент функции в произвольной точке. Определить частные производные функции f(x):
	//   grad f(x) = [ df(x)/x1,...,df(x)/dxn ]T(транспонированный вектор)
	// Начальное приближение u[0]
	x_cur = bx
	y_cur = by

	fmt.Println("## Исходная точка:")
	fmt.Printf("  x[0]: (%f, %f)\n", x_cur, y_cur)

	fmt.Println()
	fmt.Println("## Приближения:")

	// Шаг 2. Положить k = 0
	stop := false
	for k = 0; k < MAXITERATIONS && !stop ; k++ {
		// П.2. Рассчитывают vec x[k+1] = vec x[k] - lambda[k] * nabla F(vec x[k]),
		// где lambda[k]= argmin_lambda F(vec x[k] - lambda * nabla F(vec x[k]))
		//   Шаг 3. Вычислить grad f(x[k]).
		//
		//   Шаг 4. Проверить выполнение критерия окончания  ||grad f(x[k])|| < epsilon1 :
		//   а) если критерий выполнен, то x[0] = x[k], останов;
		//   б) если критерий не выполнен, то перейти к шагу 5.
		//
		//   Шаг 5. Проверить выполнение неравенства k ≥ M:
		//   а) если неравенство выполнено, то x[0] = x[k], останов;
		//   б) если нет, то перейти к шагу 6
		//
		//   Шаг 6. Вычислить величину шага lambda[k0](на начальном шаге, k = 0) из условия
		//   F(lambda[k]) = f(x[k] - lambda[k] * grad f(x[k])) -> min lambda[k]
		argmin := dichotomiaExample
		// Находим lambda[k] как минимум функции g(x[k]) на отрезке -10000,100000
		lambda = argmin(-10000, 100000, epsilon, x_cur, y_cur);
		fmt.Println("     lambdaE:", lambda)
		// Вычисляем u[k] новое прилижение
		// Шаг 7. Вычислить x[k+1] = x[k] - lambda[k] * grad f(x[k])
		x_next = x_cur - lambda * f_dx(x_cur, y_cur)
		y_next = y_cur - lambda * f_dy(x_cur, y_cur)

		fmt.Println()
		fmt.Printf("  x[%d]: (%.2f, %.2f)\n", k+1, x_next, y_next)
		fmt.Printf("  f(%.2f, %.2f) = %.2f\n", x_next, y_next, f(x_next, y_next))
		// П.3. Проверяют условие остановки:
		// Если |vec x[k+1] - vec x[k]| > epsilon
		//  или |F(vec x[k+1]) - F(vec x[k])| > epsilon
		//  или |nabla F(vec x[k+1])| > epsilon (выбирают одно из условий), 
		// то k = k + 1 и переход к П.2.
		// Иначе vec x = vec x[k+1] и останов.
		//   Шаг 8. Проверить выполнение условий ||x[k+1] - x[k]|| < epsilon2, ||f(x[k+1]) - f(x[k])|| < epsilon2:
		//   а) если оба условия выполнены при текущем значении k и k = k - 1, то расчет окончен, x[0] = x[k+1], останов;
		//   б) если хотя бы одно из условий не выполнено, то положить k = k + 1 и перейти к шагу 3
		if k > 1 {
			// Проверяем условие остановки
			if normaExample(x_next - x_cur, y_next - y_cur) < epsilon {
				// останов
				stop = true
			}
		}

		x_cur = x_next
		y_cur = y_next
	}

	fmt.Println()
	fmt.Println("## Точка минимума epsilon:", epsilon)
	fmt.Printf("   f(%.2f , %.2f) = %.2f\n", x_cur, y_cur, f(x_cur, y_cur))

// получить минимум фуркции
	return f(x_cur, y_cur)
}


//
// ТЕСТ и заодно пример
// Алгоритм моделирования одномерного стохастического процесса
func Simple_stochastic_process_model() {
	fmt.Println("simple_stochastic_process_model")

	fo,_ := os.Create("./ito.csv")
	defer fo.Close()

	SRnd64(time.Now().Unix()) // "встряхиваем" генератор

	// a(x,t) коэффициент сноса, если a(x,t) > 0, 
	// то процесс в среднем движется вверх(растёт) иначе вниз
	a := func (x, t float64) float64 { return float64(0) } // снос
	b := func (x, t float64) float64 { return x } // волатильность

	// Задаётся num = 100 значений случайного процесса x(t),
	// с равным шагом по времени: t = 0, step, 2*step,...,num*step
	// Каждый интервал между вычисляемыми значениями x разбивается 
	// на большое число lag точек, время между которыми равно dt = step/lag
	step := 0.1 // шаг по времени для табуляции х
	num := 100 // число точек табуляции
	lag := 1000 // количество дроблений шага step

	dt := float64(step) / float64(lag)
	sqrt_dt := math.Sqrt(dt) // тут хранится вычисленное заначение корня

	// исходные
	var x float64 = 1 // начальное значение x
	var t float64 = 0 // в начальный момпент времени t

	fmt.Fprintf(fo, "0;%f;%f\n", t, x)

	// Для вычислений используется два вложенных цикла.
	// Внешний по k приводит к выводу x(step), x(2step),... тогда
	// как внутренний по j, при помощи итерационной формулы, вычисляет
	// промужеточные значения случайного процесса.
	for k := 1; k <= num; k++ {
		for j := 0; j < lag; j++ {
			// новое значение x получаемое из предыдущего
			x += (a(x, t) * dt + b(x, t) * RndG()*sqrt_dt)
			// новое значение t как приращение 
			t += dt
		}

		fmt.Fprintf(fo, "%d;%f;%f\n", k, t, x)
	}
}



// ТЕСТ и заодно пример
// Двухмерный осцилятор с затуханием, имеющий скоррелированный шум <deltaWx deltaWy> = rho dt
//
// /dx = (-lambda*x - omega*y)dt + sigma*deltaWx 
//<
// \dy = (+omega*x + lambda*y)dt + sigma*deltaWy
//
func Two_dimensional_oscillator_with_attenuation() {
	// подробности на стр.252
	fmt.Println("two_dimensional_oscillator_with_attenuation")

	fo,_ := os.Create("./ito2.csv")
	defer fo.Close()

	SRnd64(time.Now().Unix()) // "встряхиваем" генератор

	// Задаётся num = 100 значений случайного процесса x(t),
	// с равным шагом по времени: t = 0, step, 2*step,...,num*step
	// Каждый интервал между вычисляемыми значениями x разбивается 
	// на большое число lag точек, время между которыми равно dt = step/lag
	step := 0.1 // шаг по времени для табуляции х
	num := 100 // число точек табуляции
	lag := 1000 // количество дроблений шага step

	w := 0.5
	lm := 0.01
	si := 1.0
	rho := 0.1

	dt := float64(step) / float64(lag)
	sqrt_dt := math.Sqrt(dt) // тут хранится вычисленное заначение корня
	rho1 := math.Sqrt(1 - rho*rho)

	// исходные
	var x_cur float64 = 1
	var y_cur float64 = 1
	var t float64 = 0
	var x_next float64
	var y_next float64
	var r1 float64
	var r2 float64

	fmt.Fprintf(fo, "0;%f;%f;%f\n", t, x_cur, y_cur)

	// Для вычислений используется два вложенных цикла.
	// Внешний по k приводит к выводу x(step), x(2step),... тогда
	// как внутренний по j, при помощи итерационной формулы, вычисляет
	// промужеточные значения случайного процесса.
	for k := 1; k <= num; k++ {
		for j := 0; j < lag; j++ {
			// если это убрать, то будут две синусойды сдвинутые по фазе :)
			r1 = RndG()
			r2 = rho * r1 + rho1*RndG()
			// новое значение x получаемое из предыдущего
			x_next = x_cur + (-lm*x_cur - w*y_cur)*dt +	si*r1*sqrt_dt
			y_next = y_cur + (-lm*y_cur + w*x_cur)*dt + si*r2*sqrt_dt
			x_cur = x_next
			y_cur = y_next
			// новое значение t как приращение 
			t += dt
		}

		fmt.Fprintf(fo, "%d;%f;%f;%f\n", k, t, x_cur, y_cur)
	}
}


// ТЕСТ и заодно пример
// пример множественной генерации
// система следующего вида:
// dxi = -xi*r*dt + deltaWi
func Multi_dimensional_oscillator_with_attenuation() {
	fmt.Println("multi_dimensional_oscillator_with_attenuation")

	fo,_ := os.Create("./ito_multi.csv")
	defer fo.Close()

	SRnd64(time.Now().Unix()) // "встряхиваем" генератор

	// Задаётся num = 100 значений случайного процесса x(t),
	// с равным шагом по времени: t = 0, step, 2*step,...,num*step
	// Каждый интервал между вычисляемыми значениями x разбивается 
	// на большое число lag точек, время между которыми равно dt = step/lag
	step := 0.1 // шаг по времени для табуляции х
	num := 300 // число точек табуляции
	lag := 1000 // количество дроблений шага step

	dt := float64(step) / float64(lag)
	sqrt_dt := math.Sqrt(dt) // тут хранится вычисленное заначение корня

	// число осциляторов
	osc_count := 200
	// исходные
	x_cur := make([]float64, osc_count)
	var t float64 = 0
	x_next := make([]float64, len(x_cur))

	// первая строка, исходная позиция
	sep := ""
	fmt.Fprintf(fo, "%s%d", sep, 0)
	sep = ";"
	for i := 0; i < len(x_cur); i++ {
		fmt.Fprintf(fo, "%s%f", sep, x_cur[i])
	}
	fmt.Fprintf(fo, "\n")

	// Для вычислений используется два вложенных цикла.
	// Внешний по k приводит к выводу x(step), x(2step),... тогда
	// как внутренний по j, при помощи итерационной формулы, вычисляет
	// промужеточные значения случайного процесса.
	for k := 1; k <= num; k++ {
		for j := 0; j < lag; j++ {
			var r float64 // расстояние до начала координат в n-мерной модели, здесь n=len(x_cur)  
			for i := 0; i < len(x_cur); i++ {
				r += x_cur[i]*x_cur[i]
			}
			r = math.Sqrt(r)
			// новое значение x получаемое из предыдущего
			for i := 0; i < len(x_cur); i++ {
				x_next[i] = x_cur[i] - x_cur[i]*r*dt + RndG()*sqrt_dt
			}

			copy(x_cur,x_next)


			// новое значение t как приращение 
			t += dt
		}
		// вывести текущую итераци в файл
		sep := ""
		fmt.Fprintf(fo, "%s%d", sep, k)
		sep = ";"
		for i := 0; i < len(x_cur); i++ {
			fmt.Fprintf(fo, "%s%f", sep, x_cur[i])
		}
		fmt.Fprintf(fo, "\n")
	}
}
*/
