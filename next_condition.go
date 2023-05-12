package main

import (
	"gonum.org/v1/gonum/mat"
)

func get_next_temp_and_load(
	ac_demand_is_ns *ScheduleData,
	brc_ot_is_n []float64,
	brm_ot_is_is_n mat.Matrix,
	brl_ot_is_is_n mat.Matrix,
	theta_lower_target_is_n []float64,
	theta_upper_target_is_n []float64,
	operation_mode_is_n []OperationMode,
	is_radiative_heating_is []bool,
	is_radiative_cooling_is []bool,
	lr_h_max_cap_is []float64,
	lr_cs_max_cap_is []float64,
	theta_natural_is_n []float64,
	n int,
) (*mat.VecDense, *mat.VecDense, *mat.VecDense) {

	var nn int
	if n < 0 {
		c := ac_demand_is_ns.Len()
		nn = n + c
	} else {
		nn = n
	}
	ac_demand_is_n := ac_demand_is_ns.Get(nn)

	roomShape := len(operation_mode_is_n)

	// 室の数
	n_room := roomShape

	// 係数　kt, W / K, [i, i], float型
	kt := brm_ot_is_is_n

	// 係数 kc, [i, i], float型
	kc := mat.NewDense(n_room, n_room, nil)
	for i := 0; i < n_room; i++ {
		kc.Set(i, i, 1)
	}

	// 係数 kr, [i, i], float型
	kr := brl_ot_is_is_n

	// 係数 k, W, [i, 1], float型
	k := brc_ot_is_n

	// 室温指定を表す係数, [i, 1], int型
	// 指定する = 0, 指定しない = 1
	// 室温を指定しない場合は、 operation_mode が STOP_CLOSE or STOP_OPEN の場合である。
	// 後で再計算する際に、負荷が機器容量を超えている場合は、最大暖房／冷房負荷で処理されることになるため、
	// 室温を指定しない場合は、この限りではない。
	nt := make([]float64, roomShape)
	for i := range nt {
		if operation_mode_is_n[i] == HEATING || operation_mode_is_n[i] == COOLING {
			nt[i] = 0 // 室温を指定する
		} else {
			nt[i] = 1 // 室温を指定しない
		}
	}

	// nt = 0 （室温を指定する） に対応する要素に、ターゲットとなるOTを代入する。
	// nt = 1 （室温を指定しない）場合は、theta_set は 0 にしなければならない。
	theta_set := mat.NewVecDense(roomShape, nil)
	for i := range nt {
		if operation_mode_is_n[i] == HEATING {
			val := theta_lower_target_is_n[i]*ac_demand_is_n[i] +
				theta_natural_is_n[i]*(1.0-ac_demand_is_n[i])
			theta_set.SetVec(i, val)
		} else if operation_mode_is_n[i] == COOLING {
			val := theta_upper_target_is_n[i]*ac_demand_is_n[i] +
				theta_natural_is_n[i]*(1.0-ac_demand_is_n[i])
			theta_set.SetVec(i, val)
		}
	}

	// 対流空調指定を表す係数, [i, 1], int型
	// 指定する = 0, 指定しない = 1
	// 対流空調を指定しない場合は、対流空調をしている場合に相当するので、
	//   operation_mode が　HEATING でかつ、 is_radiative_heating_is が false の場合か、
	//   operation_mode が COOLING でかつ、 is_radiative_cooling_is が false の場合
	// のどちらかである。
	c := make([]float64, roomShape)
	for i := range c {
		if (operation_mode_is_n[i] == HEATING && !is_radiative_heating_is[i]) ||
			(operation_mode_is_n[i] == COOLING && !is_radiative_cooling_is[i]) {
			c[i] = 1
		} else {
			c[i] = 0
		}
	}

	// c = 0 （対流空調を指定する）に対応する要素に、0.0 を代入する。
	// 対流空調を指定する場合は空調をしていないことに相当するため。ただし、後述する、最大能力で動く場合は、その値を代入することになる。
	// 対流空調を指定しない場合は、 lc_set には 0.0 を入れなければならない。
	// 結果として、ここでは、あらゆるケースで 0.0 が代入される。
	lc_set := mat.NewVecDense(roomShape, nil)

	// 放射空調指定を表す係数, [i, 1], int型
	// 指定する = 0, 指定しない = 1
	// 放射空調を指定しない場合は、放射空調をしている場合に相当するので、
	//   operation_mode が　HEATING でかつ、 is_radiative_heating_is が true の場合か、
	//   operation_mode が COOLING でかつ、 is_radiative_cooling_is が true の場合
	// のどちらかである。
	r := make([]float64, roomShape)
	for i := range r {
		if (operation_mode_is_n[i] == HEATING && is_radiative_heating_is[i]) ||
			(operation_mode_is_n[i] == COOLING && is_radiative_cooling_is[i]) {
			r[i] = 1
		} else {
			r[i] = 0
		}
	}

	// r = 0 （放射空調を指定する）に対応する要素に、0.0 を代入する。
	// 放射空調を指定する場合は空調をしていないことに相当するため。ただし、後述する、最大能力で動く場合は、その値を代入することになる。
	// 放射空調を指定しない場合は、 lr_set には 0.0 を入れなければならない。
	// 結果として、ここでは、あらゆるケースで 0.0 が代入される。
	lr_set := mat.NewVecDense(roomShape, nil)

	// theta 温度, degree C, [i, 1]
	// lc 対流空調負荷, W, [i, 1]
	// lr 放射空調負荷, W, [i, 1]
	theta, lc, lr := get_load_and_temp(kt, kc, kr, k, nt, theta_set, c, lc_set, r, lr_set)

	// 計算された放射空調負荷が最大放熱量を上回る場合は、放熱量を最大放熱量に固定して、対流空調負荷を未知数として再計算する。
	over_lr := make([]bool, roomShape)
	for i := range over_lr {
		over_lr[i] = lr.At(i, 0) > lr_h_max_cap_is[i]
	}

	for i, v := range over_lr {
		if v {
			// 対流負荷を未知数とする。
			c[i] = 1

			// 放射負荷を最大放熱量に指定する。
			r[i] = 0
			lr_set.SetVec(i, lr_h_max_cap_is[i])
		}
	}

	// 計算された放射空調負荷が最大放熱量を下回る場合は、放熱量を最大放熱量に固定して、対流空調負荷を未知数として再計算する。
	// 注意：冷房の最大放熱量は正の値で指定される。一方、計算される負荷（lr）は、冷房の場合、負の値で指定される。
	under_lr := make([]bool, roomShape)
	for i := range under_lr {
		under_lr[i] = lr.At(i, 0) < -lr_cs_max_cap_is[i]
	}

	for i, v := range under_lr {
		if v {
			// 対流負荷を未知数とする。
			c[i] = 1

			// 放射負荷を最大放熱量に指定する。
			r[i] = 0
			lr_set.SetVec(i, -lr_cs_max_cap_is[i])
		}
	}

	// 放射暖房を最大放熱量に指定して再計算する。
	theta, lc, lr = get_load_and_temp(kt, kc, kr, k, nt, theta_set, c, lc_set, r, lr_set)

	return theta, lc, lr
}

/*
   Args:
       kt: 係数 kt, W/K, [i, i], float型
       kc: 係数 kc, [i, i], float型
       kr: 係数 kr, [i, i], float型
       k: 係数 k, W, [i, 1], float型
       nt: 室温指定を表す係数, [i, 1], int型
       theta_set: 室温を指定する場合の室温, degree C, [i, 1], float型
       c: 対流空調指定を表す係数, [i, 1], int型
       lc_set: 対流空調を指定する場合の放熱量, W, [i, 1], float型
       r: 放射空調指定を表す係数, [i, 1], int型
       lr_set: 放射空調を指定する場合の放熱量, W, [i, 1], float型

   Returns:
       温度, degree C, [i, 1]
       対流空調負荷, W, [i, 1]
       放射空調負荷, W, [i, 1]

   Notes:
       各係数によって、
       kt theta = kc Lc + kr Lr + k
       が維持される。

       室温・対流空調・放射空調 を指定する場合は、指定しない = 1, 指定する = 0 とする。
       theta_set, lc_set, lr_set について、指定しない場合は必ず 0.0 とする。
*/
func get_load_and_temp(
	kt mat.Matrix,
	kc mat.Matrix,
	kr mat.Matrix,
	k []float64,
	nt []float64,
	theta_set mat.MutableVector,
	c []float64,
	lc_set mat.MutableVector,
	r []float64,
	lr_set mat.MutableVector,
) (*mat.VecDense, *mat.VecDense, *mat.VecDense) {

	n := len(nt)

	// Modify theta_set where nt == 1
	for i := 0; i < n; i++ {
		// 室温を指定しない場合 nt = 1
		// 室温を指定しない場合 theta_set は 0.0 とする。
		if nt[i] == 1 {
			theta_set.SetVec(i, 0.0)
		}

		// 対流空調負荷を指定しない場合 c = 1
		// 対流空調負荷を指定しない場合 lc_set = 0.0 とする。
		if c[i] == 1 {
			lc_set.SetVec(i, 0.0)
		}

		// 放射空調負荷を指定しない場合 r = 1
		// 放射空調負荷を指定しない場合 lr_set = 0.0 とする。
		if r[i] == 1 {
			lr_set.SetVec(i, 0.0)
		}
	}

	// Calculate x1 and x2
	ntv := mat.NewVecDense(n, nt)
	cv := mat.NewVecDense(n, c)
	rv := mat.NewVecDense(n, r)

	// kt * nt.T - kc * c.T - kr * r.T
	var x1, x1_term1, x1_term2, x1_term3 mat.Dense
	ntv_T := ntv.T()
	cv_T := cv.T()
	rv_T := rv.T()
	x1_term1.Apply(func(i, j int, v float64) float64 { return v * ntv_T.At(0, j) }, kt)
	x1_term2.Apply(func(i, j int, v float64) float64 { return v * cv_T.At(0, j) }, kc)
	x1_term3.Apply(func(i, j int, v float64) float64 { return v * rv_T.At(0, j) }, kr)
	x1.Sub(&x1_term1, &x1_term2)
	x1.Sub(&x1, &x1_term3)

	// -np.dot(kt, theta_set) + np.dot(kc, lc_set) + np.dot(kr, lr_set) + k
	var x2, x2_term1, x2_term2, x2_term3 mat.VecDense
	x2_term1.MulVec(kt, theta_set)
	x2_term2.MulVec(kc, lc_set)
	x2_term3.MulVec(kr, lr_set)
	x2.SubVec(&x2_term2, &x2_term1)
	x2.AddVec(&x2, &x2_term3)
	x2.AddVec(&x2, mat.NewVecDense(len(k), k))

	var v mat.VecDense
	v.SolveVec(&x1, &x2)

	// 求めるべき数値
	// nt, c, r それぞれ、1の場合（値を指定しない場合）は、vで表される値が入る。
	// 反対に、 0 の場合（値を指定する場合）、は、それぞれ、theta_set, lc_set, lr_set の値が入る。
	theta_rq, lc_rq, lr_rq := make([]float64, n), make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		theta_rq[i] = v.AtVec(i)*nt[i] + theta_set.AtVec(i)*(1-nt[i])
		lc_rq[i] = v.AtVec(i)*c[i] + lc_set.AtVec(i)*(1-c[i])
		lr_rq[i] = v.AtVec(i)*r[i] + lr_set.AtVec(i)*(1-r[i])
	}

	return mat.NewVecDense(n, theta_rq), mat.NewVecDense(n, lc_rq), mat.NewVecDense(n, lr_rq)
}
