package main

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

/*
傾斜面の方位角・傾斜角に応じて傾斜面の日射量を計算する。
    Args:
        w: Weather Class
        drct_j: Direction Class
    Returns:
        (1) ステップnにおける境界jの傾斜面に入射する日射量のうち直達成分, W/m2, [N+1]
        (2) ステップnにおける境界jの傾斜面に入射する日射量のうち天空成分, W/m2, [N+1]
        (3) ステップnにおける境界jの傾斜面に入射する日射量のうち地盤反射成分, W/m2, [N+1]
        (4) ステップnにおける境界jの傾斜面の夜間放射量, W/m2, [N+1]
*/
func get_i_is_j_ns(
	w *Weather,
	direction Direction,
) (i_srf_dn_j_ns, i_srf_sky_j_ns, i_srf_ref_j_ns, r_srf_eff_j_ns *mat.VecDense) {

	// ステップnにおける法線面直達日射量, W/m2K [N+1]
	i_dn_ns := w.i_dn_ns_plus
	// ステップnにおける水平面天空日射量, W/m2K [N+1]
	i_sky_ns := w.i_sky_ns_plus
	// ステップnにおける夜間放射量, W/m2, [N+1]
	r_n_ns := w.r_n_ns_plus
	// ステップnにおける太陽高度, rad [N+1]
	h_sun_ns := w.h_sun_ns_plus

	// ステップnにおける太陽方位角, rad [N+1]
	a_sun_ns := w.a_sun_ns_plus

	// 境界jの傾斜面の傾斜角, rad
	beta_w_j := direction.beta_w_j()

	// ステップnにおける境界jに入射する太陽の入射角, deg, [N+1]
	phi_j_ns := get_phi_j_ns(h_sun_ns, a_sun_ns, direction)

	// ステップ n における水平面全天日射量, W/m2, [n]
	i_hrz_ns := _get_i_hrz_ns(i_dn_ns, i_sky_ns, h_sun_ns)

	// 境界jの傾斜面の天空に対する形態係数, -
	f_sky_j := _get_f_sky_j(beta_w_j)

	// 境界 j の地面に対する傾斜面の形態係数, -
	f_gnd_j := _get_f_gnd_j(f_sky_j)

	// ステップ n における境界 j の傾斜面の夜間放射量, W/m2, [N+1]
	r_srf_eff_j_ns = _get_r_srf_eff_j_ns(r_n_ns, f_sky_j)

	// ステップnにおける境界jに入射する日射量の地盤反射成分, W/m2, [N+1]
	i_srf_ref_j_ns = _get_i_srf_ref_j_ns(f_gnd_j, i_hrz_ns)

	// ステップnにおける境界jに入射する日射量の天空成分, W/m2, [N+1]
	i_srf_sky_j_ns = _get_i_srf_sky_j_ns(i_sky_ns, f_sky_j)

	// ステップnにおける境界jに入射する日射量の直達成分, W/m2, [N+1]
	i_srf_dn_j_ns = _get_i_srf_dn_j_ns(i_dn_ns, phi_j_ns)

	return i_srf_dn_j_ns, i_srf_sky_j_ns, i_srf_ref_j_ns, r_srf_eff_j_ns
}

/*
   傾斜面に入射する日射量の直達成分を計算する。

   Args:
       i_dn_ns: ステップ n における法線面直達日射量, W/m2, [n]
       phi_j_ns: ステップ n における境界 j の傾斜面に入射する日射の入射角, rad, [n]

   Returns:
       ステップ　n　における境界 j の傾斜面に入射する日射量の直達成分, W/m2, [n]

   Notes:
       式(1)
*/
func _get_i_srf_dn_j_ns(i_dn_ns, phi_j_ns *mat.VecDense) *mat.VecDense {
	cos := mat.NewVecDense(phi_j_ns.Len(), nil)
	for i := 0; i < phi_j_ns.Len(); i++ {
		cos.SetVec(i, math.Cos(phi_j_ns.AtVec(i)))
	}

	var i_srf_dn_j_ns mat.VecDense
	i_srf_dn_j_ns.MulElemVec(i_dn_ns, cos)

	return &i_srf_dn_j_ns
}

/*
   傾斜面に入射する日射量の天空成分を計算する。

   Args:
       i_sky_ns: ステップ　n における水平面天空日射量, W/m2, [n]
       f_sky_j: 境界 j の天空に対する傾斜面の形態係数

   Returns:
       ステップ n における境界 j の傾斜面に入射する日射量の天空成分, W/m2, [n]

   Notes:
       式(2)
*/
func _get_i_srf_sky_j_ns(i_sky_ns *mat.VecDense, f_sky_j float64) *mat.VecDense {
	var i_srf_sky_j_ns mat.VecDense
	i_srf_sky_j_ns.ScaleVec(f_sky_j, i_sky_ns)
	return &i_srf_sky_j_ns
}

/*
   傾斜面の日射量のうち地盤反射成分を求める。

   Args:
       f_gnd_j: 境界 j の地面に対する傾斜面の形態係数, -
       i_hrz_ns: ステップ n における水平面全天日射量, W/m2, [n]

   Returns:
       ステップ n における境界 j の傾斜面に入射する日射量の地盤反射成分, W/m2, [n]

   Notes:
       式(3)
*/
func _get_i_srf_ref_j_ns(f_gnd_j float64, i_hrz_ns *mat.VecDense) *mat.VecDense {
	// 地面の日射反射率
	const rho_gnd = 0.1

	var i_srf_ref_j_ns mat.VecDense
	i_srf_ref_j_ns.ScaleVec(f_gnd_j*rho_gnd, i_hrz_ns)

	return &i_srf_ref_j_ns
}

/*
   傾斜面の方位角・傾斜角に応じて傾斜面の夜間放射量を計算する。

   Args:
       r_eff_ns: ステップ n における水平面の夜間放射量, W/m2, [n]
       f_sky_j: 境界 j の天空に対する傾斜面の形態係数, -

   Returns:
       ステップ n における境界 j の傾斜面の夜間放射量, W/m2, [n]

   Notes:
       式(4)
*/
func _get_r_srf_eff_j_ns(r_eff_ns *mat.VecDense, f_sky_j float64) *mat.VecDense {
	var r_srf_eff_j_ns mat.VecDense
	r_srf_eff_j_ns.ScaleVec(f_sky_j, r_eff_ns)
	return &r_srf_eff_j_ns
}

/*
   地面に対する傾斜面の形態係数を計算する。

   Args:
       f_sky_j: 境界 j の天空に対する傾斜面の形態係数, -

   Returns:
       境界 j の地面に対する傾斜面の形態係数, -

   Notes:
       式(5)
*/
func _get_f_gnd_j(f_sky_j float64) float64 {
	f_gnd_j := 1.0 - f_sky_j
	return f_gnd_j
}

/*
   傾斜面の天空に対する形態係数を計算する。

   Args:
       beta_w_j: 境界 j の傾斜面の傾斜角, rad

   Returns:
       境界jの傾斜面の天空に対する形態係数

   Notes:
       式(6)
       境界jの傾斜面の傾斜角 は水平面を0とし、垂直面をπ/2とし、オーバーハング床等における下に向いた面はπとし、値は0～πの範囲をとる。
*/
func _get_f_sky_j(beta_w_j float64) float64 {
	if beta_w_j < 0 {
		panic("傾斜面の傾斜角が0より小さい値となっています")
	}

	if beta_w_j > math.Pi {
		panic("傾斜角の傾斜面がπより大きい値となっています")
	}

	return (1.0 + math.Cos(beta_w_j)) / 2.0
}

/*
   水平面全天日射量を計算する。

   Args:
       i_dn_ns: ステップ n における法線面直達日射量, W/m2, [n]
       i_sky_ns: ステップ n における水平面天空日射量, W/m2, [n]
       h_sun_ns: ステップ n における太陽高度, rad, [n]

   Returns:
       ステップ n における水平面全天日射量, W/m2, [n]

   Notes:
       式(7)
*/
func _get_i_hrz_ns(i_dn_ns, i_sky_ns, h_sun_ns *mat.VecDense) *mat.VecDense {
	sin_h_sun_ns := mat.NewVecDense(i_dn_ns.Len(), nil)
	for i := 0; i < i_dn_ns.Len(); i++ {
		h_sun := math.Max(h_sun_ns.AtVec(i), 0)
		sin_h_sun_ns.SetVec(i, math.Sin(h_sun))
	}

	var i_hsr_ns mat.VecDense
	i_hsr_ns.MulElemVec(sin_h_sun_ns, i_dn_ns)
	i_hsr_ns.AddVec(&i_hsr_ns, i_sky_ns)

	return &i_hsr_ns
}

/*
   傾斜面に入射する太陽の入射角を計算する。

   Args:
        h_sun_ns: ステップnにおける太陽高度, rad, [N+1]
        a_sun_ns: ステップnにおける太陽方位角, rad, [N+1]
        drct_j: Direction Class

   Returns:
       ステップ n の境界 j における傾斜面に入射する太陽の入射角, rad, [n]

   Notes:
       式(8), 式(9)
*/
func get_phi_j_ns(h_sun_ns, a_sun_ns *mat.VecDense, drct_j Direction) *mat.VecDense {
	theta_aoi_j_ns := mat.NewVecDense(h_sun_ns.Len(), nil)

	if drct_j == DirectionTop {
		// 方位が上面（beta_w_j=0）の場合は、厳密には方位角（alpha_w_j）は定義できないため、
		// 条件分岐により式を分ける。
		for i := 0; i < h_sun_ns.Len(); i++ {
			var cos_phi_j_ns float64
			h_sun := h_sun_ns.At(i, 0)
			cos_phi_j_ns = math.Max(math.Sin(h_sun), 0)
			// ステップnにおける境界jに入射する日射の入射角, rad, [N+1]
			theta_aoi_j_ns.SetVec(i, math.Acos(cos_phi_j_ns))
		}
	} else if drct_j == DirectionBottom {
		// 方位が下面（beta_w_j=0）の場合は、厳密には方位角（alpha_w_j）は定義できないため、
		// 条件分岐により式を分ける。
		for i := 0; i < h_sun_ns.Len(); i++ {
			// ステップnにおける境界jに入射する日射の入射角, rad, [N+1]
			const acos_0 = math.Pi / 2
			theta_aoi_j_ns.SetVec(i, acos_0)
		}
	} else {
		// 境界 j における傾斜面の方位角, rad
		alpha_w_j := drct_j.alpha_w_j()

		// 境界 j における傾斜面の傾斜角, rad
		beta_w_j := drct_j.beta_w_j()
		cos_beta := math.Cos(beta_w_j)
		sin_beta := math.Sin(beta_w_j)

		for i := 0; i < h_sun_ns.Len(); i++ {
			var cos_phi_j_ns float64
			h_sun := h_sun_ns.At(i, 0)

			// ステップ n の境界 j における傾斜面に入射する太陽の入射角の余弦, -, [n]
			// cos(h_sun_ns) == 0.0 の場合は太陽が天頂にある時であり、太陽の方位角が定義されない。
			// その場合、cos(h_sun_ns)がゼロとなり、下式の第2項・第3項がゼロになる。
			// これを回避するために場合分けを行っている。
			// 余弦がマイナス（入射角が90°～270°）の場合は傾斜面の裏面に太陽が位置していることになるため、値をゼロにする。
			// （法線面直達日射量にこの値をかけるため、結果的に日射があたらないという計算になる。）
			cos_h_sun := math.Cos(h_sun)
			sin_h_sun := math.Sin(h_sun)
			if cos_h_sun == 0.0 {
				cos_phi_j_ns = math.Max(sin_h_sun*cos_beta, 0)
			} else {
				a_sun := a_sun_ns.AtVec(i)
				cos_phi_j_ns = math.Max(sin_h_sun*cos_beta+
					cos_h_sun*math.Sin(a_sun)*sin_beta*math.Sin(alpha_w_j)+
					cos_h_sun*math.Cos(a_sun)*sin_beta*math.Cos(alpha_w_j), 0)
			}

			// ステップnにおける境界jに入射する日射の入射角, rad, [N+1]
			theta_aoi_j_ns.SetVec(i, math.Acos(cos_phi_j_ns))
		}
	}

	return theta_aoi_j_ns
}
