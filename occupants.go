package main

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
	1人あたりの人体発湿を計算する。

    Args:
        theta_r_is_n: ステップnの室iにおける室温, degree C, [i, 1]

    Returns:
        ステップnの室iにおける1人あたりの人体発湿, kg/s, [i, 1]
*/
func get_x_hum_psn_is_n(theta_r_is_n []float64, q_hum_psn_is_n []float64) []float64 {
	x_hum_psn_is_n := make([]float64, len(theta_r_is_n))
	for i := 0; i < len(theta_r_is_n); i++ {
		x_hum_psn_is_n[i] = (119.0 - q_hum_psn_is_n[i]) / l_wtr
	}

	return x_hum_psn_is_n
}

/*
	1人あたりの人体発湿を計算する。

    Args:
        theta_r_is_n: ステップnの室iにおける室温, degree C, [i, 1]

    Returns:
        ステップnの室iにおける1人あたりの人体発熱, W, [i, 1]
*/
func get_q_hum_psn_is_n(theta_r_is_n mat.Vector) mat.Vector {
	q_hum_psn_is_n := mat.NewVecDense(theta_r_is_n.Len(), nil)
	for i := 0; i < theta_r_is_n.Len(); i++ {
		q_hum_psn_is_n.SetVec(i, math.Max(63.0-4.0*(theta_r_is_n.AtVec(i)-24.0), 119.0))
	}
	return q_hum_psn_is_n
}

func get_q_hum_psn_is_n_slice(theta_r_is_n []float64) []float64 {
	q_hum_psn_is_n := make([]float64, len(theta_r_is_n))
	for i := 0; i < len(theta_r_is_n); i++ {
		q_hum_psn_is_n[i] = math.Max(63.0-4.0*(theta_r_is_n[i]-24.0), 119.0)
	}
	return q_hum_psn_is_n
}

const clo_heavy = 1.1

const clo_middle = 0.7

const clo_light = 0.3
