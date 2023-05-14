package main

import (
	"gonum.org/v1/gonum/mat"
)

/*
   境界jが接する室の在室者に対する境界jの形態係数を取得する。
   Args:
       p_is_js: 室と境界との関係を表すベクトル
       a_s_js: 境界 j の面積, m2, [j, 1]
       is_floor_js: 境界 j が床かどうか, [j, 1]
   Returns:
       境界jが接する室の在室者に対する境界jの形態係数
*/
func get_f_mrt_hum_js(p_is_js *mat.Dense, a_s_js mat.Vector, is_floor_js []bool) *mat.Dense {
	r, c := p_is_js.Dims()

	a_s_js_floor := make([]float64, a_s_js.Len())
	a_s_js_not_floor := make([]float64, a_s_js.Len())
	for i := 0; i < a_s_js.Len(); i++ {
		if is_floor_js[i] {
			// boundaries area which is floor, m2, [j, 1]
			a_s_js_floor[i] = a_s_js.AtVec(i)
		} else {
			// boundaries area which is not floor, m2, [j, 1]
			a_s_js_not_floor[i] = a_s_js.AtVec(i)
		}
	}

	var temp1, temp2 mat.VecDense
	temp1.MulVec(p_is_js, mat.NewVecDense(len(a_s_js_floor), a_s_js_floor))
	temp2.MulVec(p_is_js, mat.NewVecDense(len(a_s_js_not_floor), a_s_js_not_floor))

	_temp1 := temp1.RawVector().Data[:r]
	_temp2 := temp2.RawVector().Data[:r]

	f_mrt_hum_is_js := make([]float64, r*c)
	off := 0
	for i := 0; i < r; i++ {
		_p_is_js := p_is_js.RawRowView(i)
		for j := 0; j < c; j++ {
			v1 := _p_is_js[j] * a_s_js_floor[j] / _temp1[i] * 0.45
			v2 := _p_is_js[j] * a_s_js_not_floor[j] / _temp2[i] * 0.55
			f_mrt_hum_is_js[off] = v1 + v2
			off++
		}
	}

	return mat.NewDense(r, c, f_mrt_hum_is_js)
}
