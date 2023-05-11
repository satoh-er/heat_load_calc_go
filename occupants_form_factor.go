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
func get_f_mrt_hum_js(p_is_js mat.Matrix, a_s_js mat.Vector, is_floor_js []bool) *mat.Dense {
	r, c := p_is_js.Dims()

	a_s_js_floor := mat.NewVecDense(a_s_js.Len(), nil)
	a_s_js_not_floor := mat.NewVecDense(a_s_js.Len(), nil)
	for i := 0; i < a_s_js.Len(); i++ {
		if is_floor_js[i] {
			// boundaries area which is floor, m2, [j, 1]
			a_s_js_floor.SetVec(i, a_s_js.AtVec(i))
		} else {
			// boundaries area which is not floor, m2, [j, 1]
			a_s_js_not_floor.SetVec(i, a_s_js.AtVec(i))
		}
	}

	var temp1, temp2 mat.VecDense
	temp1.MulVec(p_is_js, a_s_js_floor)
	temp2.MulVec(p_is_js, a_s_js_not_floor)
	f_mrt_hum_is_js := mat.NewDense(r, c, nil)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			v1 := p_is_js.At(i, j) * a_s_js_floor.AtVec(j) / temp1.AtVec(i) * 0.45
			v2 := p_is_js.At(i, j) * a_s_js_not_floor.AtVec(j) / temp2.AtVec(i) * 0.55
			f_mrt_hum_is_js.Set(i, j, v1+v2)
		}
	}

	return f_mrt_hum_is_js
}
