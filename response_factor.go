package main

// 応答係数の初項、指数項別応答係数、公比の計算

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

const nRoot = 12

type ResponseFactor struct {
	_rft0 float64       // 貫流応答係数の初項
	_rfa0 float64       // 吸熱応答係数の初項
	_rft1 *mat.VecDense // 貫流応答係数
	_rfa1 *mat.VecDense // 吸熱応答係数
	_row  *mat.VecDense // 公比
}

func NewResponseFactor(
	rft0 float64,
	rfa0 float64,
	rft1 *mat.VecDense,
	rfa1 *mat.VecDense,
	row *mat.VecDense,
) *ResponseFactor {
	return &ResponseFactor{
		_rft0: rft0,
		_rfa0: rfa0,
		_rft1: rft1,
		_rfa1: rfa1,
		_row:  row,
	}
}

// 貫流応答係数の初項
func (rf *ResponseFactor) rft0() float64 {
	return rf._rft0
}

// 吸熱応答係数の初項
func (rf *ResponseFactor) rfa0() float64 {
	return rf._rfa0
}

// 貫流応答係数
func (rf *ResponseFactor) rft1() *mat.VecDense {
	return rf._rft1
}

// 吸熱応答係数
func (rf *ResponseFactor) rfa1() *mat.VecDense {
	return rf._rfa1
}

// 公比
func (rf *ResponseFactor) row() *mat.VecDense {
	return rf._row
}

/*
	Args:
		u_w: 熱貫流率, W/m2K
		r_i: 室内側（総合）熱伝達抵抗, m2K/W

	Returns:
		応答係数
*/
func create_for_steady(u_w float64, r_i float64) *ResponseFactor {
	// 開口部の室内表面から屋外までの熱コンダクタンス, W/m2K
	u_so := 1.0 / (1.0/u_w - r_i)

	return NewResponseFactor(
		1.0,
		1.0/u_so,
		mat.NewVecDense(nRoot, nil),
		mat.NewVecDense(nRoot, nil),
		mat.NewVecDense(nRoot, nil),
	)
}

/*
	応答係数を作成する（地盤以外に用いる）

	裏面に、熱抵抗をもち、熱容量は 0.0 の層を追加する。
	Args:
		cs: 単位面積あたりの熱容量, kJ/m2K, [layer数]
		rs: 熱抵抗, m2K/W, [layer数]
		r_o: 室外側熱伝達抵抗, m2K/W

	Returns:
		応答係数
*/
func create_for_unsteady_not_ground(cs []float64, rs []float64, r_o float64) *ResponseFactor {
	// 裏面に熱容量 0.0 、熱抵抗 r_o の層を加える。
	cs = append(cs[:], 0.0)
	rs = append(rs[:], r_o)

	// 単位変換 kJ/m2K -> J/m2K
	for i := range cs {
		cs[i] *= 1000.0
	}

	// 応答係数
	rft0, rfa0, rft1, rfa1, row := calc_response_factor_non_residential(cs, rs)

	return NewResponseFactor(rft0, rfa0, rft1, rfa1, row)
}

/*
	応答係数を作成する（地盤用）

	地盤層を追加する。
	Args:
		cs: 単位面積あたりの熱容量, kJ/m2K, [layer数]
		rs: 熱抵抗, m2K/W, [layer数]

	Returns:
		応答係数
*/
func create_for_unsteady_ground(cs []float64, rs []float64) *ResponseFactor {
	// 裏面に地盤の層を加える。
	cs = append(cs[:], 3300.0*3.0)
	rs = append(rs[:], 3.0/1.0)

	// 単位変換 kJ/m2K -> J/m2K
	for i := range cs {
		cs[i] *= 1000.0
	}

	// 応答係数
	rft0, rfa0, rft1, rfa1, row := calc_response_factor(true, cs, rs)

	// 貫流応答係数の上書
	// 土壌の計算は吸熱応答のみで計算するため、畳み込み積分に必要な指数項別応答係数はすべて０にする
	// 貫流応答の初項は年平均気温掛かる係数であることから１とし、計算された貫流応答係数をすべて上書きする
	rft0 = 1.0
	rft1 = mat.NewVecDense(nRoot, nil)

	return NewResponseFactor(rft0, rfa0, rft1, rfa1, row)
}

// ラプラス変数の設定
// alp: 固定根
func get_laps(alp []float64) []float64 {
	// 与えるラプラス変数の個数
	n := len(alp) * 2

	// ラプラス変数の配列
	laps := make([]float64, n)

	for i := 1; i <= n; i++ {
		if i%2 == 0 {
			// 偶数番目はαをそのまま入力
			laps[i-1] = alp[int((i-1)/2)]
		} else if i == 1 {
			// 最初はα1/√(α2/α1）とする
			laps[i-1] = alp[0] / math.Sqrt(alp[1]/alp[0])
		} else {
			// それ以外は等比数列で補間
			lngL := int(math.Ceil(float64(i-1) / 2))
			laps[i-1] = alp[lngL] / math.Sqrt(alp[lngL]/alp[lngL-1])
		}
	}

	return laps
}

/*
   固定根を取得する。

   Args:
       is_ground: 熱容量を持つ外皮が地盤かどうか。（地盤でない場合は、一般部位または間仕切り）

   Returns:
       固定根
*/
func get_alpha_m(is_ground bool) []float64 {
	if is_ground {
		// 地盤の場合
		return []float64{
			1.05699306612549e-08,
			3.27447457677204e-08,
			1.01440436059147e-07,
			3.14253839100315e-07,
			9.73531652917036e-07,
			3.01591822058485e-06,
			9.34305801562961e-06,
			2.89439987091205e-05,
			8.96660450863221e-05,
			2.77777777777778e-04,
		}
	} else {
		// 地盤以外の場合
		return []float64{
			2.0000e-06,
			7.9000e-06,
			3.1205e-05,
			1.2325e-04,
			4.8687e-04,
			1.9231e-03,
			7.5964e-03,
			3.0006e-02,
		}
	}
}

func _sum(a []float64) float64 {
	var sum float64
	for _, v := range a {
		sum += v
	}
	return sum
}

// 壁体の単位応答の計算
// Args:
// 		C_i_k_p: 壁体の熱容量, kJ/m2K, [layer数]
// 		R_i_k_p: 壁体の熱抵抗, m2K/W, [layer数]
// 		laps : ラプラス変数
// 		alp : 固定根
// 		M : 積分の上限
func get_step_reps_of_wall(C_i_k_p, R_i_k_p, laps, alp []float64, M int) (float64, float64, []float64, []float64) {
	// 吸熱、貫流の各伝達関数ベクトルの初期化
	nlaps := len(laps)
	matGA := mat.NewVecDense(nlaps, nil)
	matGT := mat.NewVecDense(nlaps, nil)

	// 吸熱、貫流の各伝達関数ベクトルの初期化
	dblAT0 := 1.0
	dblAA0 := _sum(R_i_k_p)

	// GA(0), GT(0)
	dblGA0 := dblAA0
	dblGT0 := dblAT0

	// 壁体の熱容量が0（定常）の場合
	// 定常部位であっても、そのまま処理を継続する（計算上は問題ないため）
	// if abs(dblCtotal) < 0.001:
	//    pass //　暫定処理（VBAではここで処理を抜ける）

	// 吸熱、貫流の各伝達関数ベクトルの作成
	for lngI := 0; lngI < len(laps); lngI++ {
		// 伝達関数の計算
		GA, GT := calc_transfer_function(C_i_k_p, R_i_k_p, laps[lngI])

		// 吸熱、貫流の各伝達関数ベクトルの作成
		matGA.SetVec(lngI, GA-dblGA0)
		matGT.SetVec(lngI, GT-dblGT0)
	}

	// 伝達関数の係数を求めるための左辺行列を作成
	nroot := len(alp)
	matF := mat.NewDense(nlaps, nroot, nil)
	for lngI, lap := range laps {
		for lngJ, root := range alp {
			matF.Set(lngI, lngJ, lap/(lap+root))
		}
	}

	// 最小二乗法のための係数行列を作成
	var matU mat.Dense
	matU.Mul(matF.T(), matF)

	// 最小二乗法のための定数項行列を作成
	var matFGA, matFGT mat.Dense
	matFGA.Apply(func(i, j int, v float64) float64 { return v * matGA.AtVec(i) }, matF)
	matFGT.Apply(func(i, j int, v float64) float64 { return v * matGT.AtVec(i) }, matF)
	CA := make([]float64, nroot)
	CT := make([]float64, nroot)
	for lngI := 0; lngI < nroot; lngI++ {
		for lngJ := 0; lngJ < nlaps; lngJ++ {
			CA[lngI] += matFGA.At(lngJ, lngI)
			CT[lngI] += matFGT.At(lngJ, lngI)
		}
	}
	matCA := mat.NewVecDense(nroot, CA)
	matCT := mat.NewVecDense(nroot, CT)

	// 伝達関数の係数を計算
	var matAA, matAT mat.VecDense
	matAA.SolveVec(&matU, matCA)
	matAT.SolveVec(&matU, matCT)

	// 伝達関数の係数を一次元配列に変換
	dblAT := matToSlice(&matAT)
	dblAA := matToSlice(&matAA)

	return dblAT0, dblAA0, dblAT, dblAA
}

/*
	伝達関数の計算
    Args:
        C_i_k_p: 層の熱容量[J/m2K]
        R_i_k_p: 層の熱抵抗[m2K/W]
        laps: ラプラス変数[1/s]

    Returns:
        貫流伝達関数[K]
*/
func calc_transfer_function(C_i_k_p []float64, R_i_k_p []float64, laps float64) (float64, float64) {
	// 四端子行列の初期化
	matFt := mat.NewDense(2, 2, []float64{
		1.0, 0.0,
		0.0, 1.0,
	})
	var matFi *mat.Dense
	for lngK, R_k := range R_i_k_p {
		C_k := C_i_k_p[lngK]

		// ---- 四端子基本行列 matFi ----
		if math.Abs(C_k) < 0.001 {
			// 定常部位（空気層等）の場合
			matFi = mat.NewDense(2, 2, []float64{
				1.0, R_k,
				0.0, 1.0,
			})
		} else {
			// 非定常部位の場合
			dblTemp := math.Sqrt(R_k * C_k * laps)
			dblCosh := math.Cosh(dblTemp)
			dblSinh := math.Sinh(dblTemp)

			matFi = mat.NewDense(2, 2, []float64{
				dblCosh, R_k / dblTemp * dblSinh,
				dblTemp / R_k * dblSinh, dblCosh,
			})
		}

		// ---- 四端子行列 matFt ----
		matFt.Mul(matFt, matFi)
	}

	// 吸熱、貫流の各伝達関数ベクトルの作成
	GA := matFt.At(0, 1) / matFt.At(1, 1)
	GT := 1.0 / matFt.At(1, 1)

	return GA, GT
}

// 壁体の単位応答の計算（非住宅向け重み付き最小二乗法適用）
// Args:
//	C_i_k_p: 熱容量の配列
//	R_i_k_p: 熱抵抗の配列
//	laps: 時刻の配列
//	alp: 温度時定数の配列
//	M: 応答係数で作成する項数
func get_step_reps_of_wall_weighted(C_i_k_p, R_i_k_p []float64, laps []float64, alp []float64, M int) (float64, float64, []float64, []float64, []float64, []float64) {
	// 四端子基本行列の初期化
	matFi := make([]*mat.Dense, len(C_i_k_p))

	// 吸熱、貫流の各伝達関数ベクトルの初期化
	nlaps := len(laps)
	matGA := mat.NewVecDense(nlaps, nil)
	matGT := mat.NewVecDense(nlaps, nil)

	// 単位貫流応答、単位吸熱応答の初期化
	dblAT0 := 1.0
	dblAA0 := _sum(R_i_k_p)

	// GA(0), GT(0)
	dblGA0 := dblAA0
	dblGT0 := dblAT0

	// 吸熱、貫流の各伝達関数ベクトルの作成
	matFt := mat.NewDense(2, 2, nil)
	for lngI := 0; lngI < len(laps); lngI++ {
		for lngK, R_k := range R_i_k_p {
			C_k := C_i_k_p[lngK]

			// ---- 四端子基本行列 matFi ----
			var matFiElem *mat.Dense
			if math.Abs(C_k) < 0.001 {
				// 定常部位（空気層等）の場合
				matFiElem = mat.NewDense(2, 2, []float64{
					1.0, R_k,
					0.0, 1.0,
				})
			} else {
				// 非定常部位の場合
				dblTemp := math.Sqrt(R_k * C_k * laps[lngI])
				dblCosh := math.Cosh(dblTemp)
				dblSinh := math.Sinh(dblTemp)

				matFiElem = mat.NewDense(2, 2, []float64{
					dblCosh, R_k / dblTemp * dblSinh,
					dblTemp / R_k * dblSinh, dblCosh,
				})
			}

			matFi[lngK] = matFiElem

			// ---- 四端子行列 matFt ----
			if lngK == 0 {
				// 室内側1層目の場合は、四端子行列に四端子基本行列をコピーする
				matFt.CloneFrom(matFiElem)
			} else {
				// 室内側2層目以降は、四端子基本行列を乗算
				matFt.Mul(matFt, matFiElem)
			}
		}

		// 吸熱、貫流の各伝達関数ベクトルの作成
		matGA.SetVec(lngI, matFt.At(0, 1)/matFt.At(1, 1)-dblGA0)
		matGT.SetVec(lngI, 1.0/matFt.At(1, 1)-dblGT0)
	}

	// 伝達関数の係数を求めるための左辺行列を作成
	nroot := len(alp)
	matF := mat.NewDense(nlaps, nroot, nil)
	for lngI, lap := range laps {
		for lngJ, root := range alp {
			matF.Set(lngI, lngJ, lap/(lap+root))
		}
	}

	// 最小二乗法のための係数行列を作成
	matU := mat.NewDense(nroot, nroot, nil)
	for lngK := 0; lngK < nroot; lngK++ {
		for lngJ := 0; lngJ < nroot; lngJ++ {
			for lngI := 0; lngI < nlaps; lngI++ {
				laps_2 := laps[lngI] * laps[lngI]
				matU.Set(lngK, lngJ, matU.At(lngK, lngJ)+laps_2*matF.At(lngI, lngK)*matF.At(lngI, lngJ))
			}
		}
	}

	// 最小二乗法のための定数項行列を作成
	matCA := mat.NewVecDense(nroot, nil)
	matCT := mat.NewVecDense(nroot, nil)
	for lngK := 0; lngK < nroot; lngK++ {
		for lngI := 0; lngI < nlaps; lngI++ {
			laps_2_F := laps[lngI] * laps[lngI] * matF.At(lngI, lngK)
			matCA.SetVec(lngK, matCA.At(lngK, 0)+laps_2_F*matGA.At(lngI, 0))
			matCT.SetVec(lngK, matCT.At(lngK, 0)+laps_2_F*matGT.At(lngI, 0))
		}
	}

	var matU_inv mat.Dense
	matU_inv.Inverse(matU)

	// 伝達関数の係数を計算
	var matAA, matAT mat.VecDense
	// err1 := matAA.SolveVec(matU, matCA)
	// err2 := matAT.SolveVec(matU, matCT)
	matAA.MulVec(&matU_inv, matCA)
	matAT.MulVec(&matU_inv, matCT)

	// if err1 != nil || err2 != nil {
	// 	panic("Error during matrix solving")
	// }

	// 伝達関数の係数を一次元配列に変換
	dblAT := matToSlice(&matAT)
	dblAA := matToSlice(&matAA)

	dblATstep := make([]float64, M)
	dblAAstep := make([]float64, M)
	for lngJ := 0; lngJ < M; lngJ++ {
		dblATstep[lngJ] = dblAT0
		dblAAstep[lngJ] = dblAA0
		for lngK, root := range alp {
			eterm := math.Exp(-root * float64(lngJ) * 900)
			dblATstep[lngJ] += dblAT[lngK] * eterm
			dblAAstep[lngJ] += dblAA[lngK] * eterm
		}
	}

	return dblAT0, dblAA0, dblAT, dblAA, dblATstep, dblAAstep
}

func matToSlice(matrix *mat.VecDense) []float64 {
	rows, cols := matrix.Dims()
	slice := make([]float64, rows*cols)
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			slice[i*cols+j] = matrix.At(i, j)
		}
	}
	return slice
}

// 二等辺三角波励振の応答係数、指数項別応答係数、公比の計算
func get_RFTRI(alp []float64, AT0 float64, AA0 float64, AT []float64, AA []float64, M int) ([]float64, []float64, []float64, []float64, []float64) {
	// 二等辺三角波励振の応答係数の配列を初期化
	dblRFT := make([]float64, M)
	dblRFA := make([]float64, M)

	// 二等辺三角波励振の応答係数の初項を計算
	dblTemp := make([]float64, len(alp))
	for i, val := range alp {
		dblTemp[i] = val * 900
	}

	// 二等辺三角波励振の応答係数の二項目以降を計算
	dblE1 := make([]float64, len(alp))
	for i, val := range dblTemp {
		dblE1[i] = (1.0 - math.Exp(-val)) / val
	}
	dblRFT[0] = AT0 + _sum(multiply(dblE1, AT))
	dblRFA[0] = AA0 + _sum(multiply(dblE1, AA))

	// 二等辺三角波励振の応答係数の二項目以降を計算
	for lngJ := 1; lngJ < M; lngJ++ {
		for i, val := range dblTemp {
			temp := (1.0 - math.Exp(-val))
			dblE1[i] = temp * temp * math.Exp(-(float64(lngJ)-1.0)*val) / val
		}
		dblRFT[lngJ] = -_sum(multiply(dblE1, AT))
		dblRFA[lngJ] = -_sum(multiply(dblE1, AA))
	}

	// 指数項別応答係数、公比を計算
	for i, val := range dblTemp {
		dblE1[i] = (1.0 - math.Exp(-val)) * (1.0 - math.Exp(-val)) / val
	}
	dblRFT1 := make([]float64, len(AT))
	dblRFA1 := make([]float64, len(AA))
	for i := range AT {
		dblRFT1[i] = -AT[i] * dblE1[i]
		dblRFA1[i] = -AA[i] * dblE1[i]
	}
	dblRow := make([]float64, len(alp))
	for i, val := range dblTemp {
		dblRow[i] = math.Exp(-val)
	}

	return dblRFT, dblRFA, dblRFT1, dblRFA1, dblRow
}

func multiply(slice1, slice2 []float64) []float64 {
	result := make([]float64, len(slice1))
	for i := range slice1 {
		result[i] = slice1[i] * slice2[i]
	}
	return result
}

// 応答係数
func calc_response_factor(is_ground bool, cs, rs []float64) (float64, float64, *mat.VecDense, *mat.VecDense, *mat.VecDense) {
	NcalTime := 50                  // 応答係数を作成する時間数[h]
	M := int(NcalTime*3600/900) + 1 //応答係数で作成する項数

	//  固定根, 一般部位の場合[8], 地盤の場合[10]
	alpha_m := get_alpha_m(is_ground)

	// ラプラス変数の設定
	laps := get_laps(alpha_m)

	// 単位応答の計算
	AT0, AA0, AT, AA := get_step_reps_of_wall(cs, rs, laps, alpha_m, M)

	// 二等辺三角波励振の応答係数、指数項別応答係数、公比の計算
	RFT, RFA, RFT1, RFA1, Row := get_RFTRI(alpha_m, AT0, AA0, AT, AA, M)

	RFT0 := RFT[0] // 貫流応答係数の初項
	RFA0 := RFA[0] // 吸熱応答係数の初項
	//Nroot := len(alpha_m) // 根の数

	RFT1_12 := make([]float64, 12)
	RFA1_12 := make([]float64, 12)
	Row_12 := make([]float64, 12)
	copy(RFT1_12, RFT1)
	copy(RFA1_12, RFA1)
	copy(Row_12, Row)

	matRFT1 := mat.NewVecDense(len(RFT1_12), RFT1_12)
	matRFA1 := mat.NewVecDense(len(RFA1_12), RFA1_12)
	matRow := mat.NewVecDense(len(Row_12), Row_12)

	return RFT0, RFA0, matRFT1, matRFA1, matRow
}

func logspace(start, stop float64, num int, base float64) []float64 {
	if num <= 0 {
		return []float64{}
	}

	step := (stop - start) / float64(num-1)
	logspaceSlice := make([]float64, num)

	for i := 0; i < num; i++ {
		logspaceSlice[i] = math.Pow(base, start+float64(i)*step)
	}

	return logspaceSlice
}

// 応答係数（非住宅用　住宅との相違は固定根と重み付き最小二乗法を使用する点）
func calc_response_factor_non_residential(C_i_k_p, R_i_k_p []float64) (float64, float64, *mat.VecDense, *mat.VecDense, *mat.VecDense) {
	NcalTime := 50                           // 応答係数を作成する時間数[h]
	M := int(float64(NcalTime*3600)/900) + 1 // 応答係数で作成する項数

	// 固定根, 初項 1/(86400*365)、終項 1/600、項数 10
	first_term := math.Log10(1.0 / (86400.0 * 365.0))
	last_term := math.Log10(1.0 / 600.0)
	alpha_m := logspace(first_term, last_term, 10, 10)

	// ラプラス変数の設定
	laps := get_laps(alpha_m)

	// 単位応答の計算
	AT0, AA0, AT, AA, _, _ := get_step_reps_of_wall_weighted(C_i_k_p, R_i_k_p, laps, alpha_m, M)

	// 二等辺三角波励振の応答係数、指数項別応答係数、公比の計算
	RFT, RFA, RFT1, RFA1, Row := get_RFTRI(alpha_m, AT0, AA0, AT, AA, M)

	RFT0 := RFT[0] // 貫流応答係数の初項
	RFA0 := RFA[0] // 吸熱応答係数の初項
	//Nroot := len(alpha_m) // 根の数

	RFT1_12 := make([]float64, 12)
	RFA1_12 := make([]float64, 12)
	Row_12 := make([]float64, 12)
	copy(RFT1_12, RFT1)
	copy(RFA1_12, RFA1)
	copy(Row_12, Row)

	matRFT1 := mat.NewVecDense(len(RFT1_12), RFT1_12)
	matRFA1 := mat.NewVecDense(len(RFA1_12), RFA1_12)
	matRow := mat.NewVecDense(len(Row_12), Row_12)

	return RFT0, RFA0, matRFT1, matRFA1, matRow
}
