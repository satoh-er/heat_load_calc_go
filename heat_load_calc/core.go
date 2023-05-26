package main

import (
	"log"
	"time"
)

/*
coreメインプログラム

    Args:
        rd: 住宅計算条件
        w: 外界気象条件
        scd: スケジュール
        itv: 時間間隔
        n_step_hourly: 計算間隔（1時間を何分割するかどうか）（デフォルトは4（15分間隔））
        n_d_main: 本計算を行う日数（デフォルトは365日（1年間））, d
        n_d_run_up: 助走計算を行う日数（デフォルトは365日（1年間））, d
        n_d_run_up_build: 助走計算のうち建物全体を解く日数（デフォルトは183日（およそ半年））, d

    Returns:
        以下のタプル
            (1) 計算結果（詳細版）をいれたDataFrame
            (2) 計算結果（簡易版）をいれたDataFrame

    Notes:
        「助走計算のうち建物全体を解く日数」は「助走計算を行う日数」で指定した値以下でないといけない。
*/
func calc(
	rd map[string]interface{},
	w *Weather,
	scd *Schedule,
	itv Interval,
	n_step_hourly int,
	n_d_main int,
	n_d_run_up int,
	n_d_run_up_build int,
) {
	log.Printf("計算開始")

	// 本計算のステップ数
	// 助走計算のステップ数
	// 助走計算のうち建物全体を解くステップ数
	n_step_main, n_step_run_up, n_step_run_up_build := get_n_step(n_step_hourly, n_d_main, n_d_run_up, n_d_run_up_build)

	// 時間間隔, s
	//deltaT := itv.get_delta_t()

	// json, csv ファイルからパラメータをロードする。
	// （ループ計算する必要の無い）事前計算を行い, クラス PreCalcParameters, PreCalcParametersGround に必要な変数を格納する。
	sqc := NewSequence(itv, rd, w, scd)

	pp := sqc.pre_calc_parameters

	gc_n := initialize_ground_conditions(sqc.bs.n_ground)

	log.Println("助走計算（土壌のみ）")
	N := scd.ac_demand_is_ns.Len()

	nn := N - n_step_run_up_build
	for n := -n_step_run_up; n < -n_step_run_up_build; n++ {
		gc_n = sqc.run_tick_ground(gc_n, n, nn)
		nn++
	}

	result := NewRecorder(n_step_main, sqc.rms.id_rm_is, sqc.bs.id_bs_js, itv)

	result.pre_recording(pp)

	// 建物を計算するにあたって初期値を与える
	c_n := initialize_conditions(sqc.rms.n_rm, sqc.bs.n_b)

	// 地盤計算の結果（項別公比法の指数項mの吸熱応答の項別成分・表面熱流）を建物の計算に引き継ぐ
	c_n = update_conditions_by_ground_conditions(sqc.bs.is_ground_js, c_n, gc_n)

	log.Println("助走計算（建物全体）")

	N_plus := N + 1
	nn = N - n_step_run_up_build
	nn_plus := N_plus - n_step_run_up_build
	for n := -n_step_run_up_build; n < 0; n++ {
		c_n, _ = sqc.run_tick(n, nn, nn_plus, c_n, result)
		nn++
	}

	log.Println("本計算")

	m := 1
	var l_cs_is_n []float64
	var l_cs_is [15]float64
	for n := 0; n < n_step_main; n++ {
		c_n, l_cs_is_n = sqc.run_tick(n, n, n, c_n, result)
		for i, x := range l_cs_is_n{
			l_cs_is[i] += x
		}
		
		if n == int(float64(n_step_main)/12*float64(m)) {
			log.Printf("%d / 12 calculated.", m)
			m++
		}
	}
	for i := 0; i < 15; i++{
		println(l_cs_is[i])
	}

	// result.post_recording(pp)

	log.Println("ログ作成")

	// dd: data detail, 15分間隔のすべてのパラメータ pd.DataFrame
	// dd_i, dd_a, err := result.export_pd()
	// if err != nil {
	// 	return nil, nil, err
	// }

	// return dd_i, dd_a, nil
}

func get_h_and_c_period(region_code int) (time.Time, time.Time, time.Time, time.Time) {

	(heating_st, heating_en, cooling_st, cooling_en) := switch region_code {
		case 1: (time.Date(1989, 9, 24, 0, 0, 0, 0, time.Local),
				time.Date(1989, 6, 7, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 7, 10, 0, 0, 0, 0, time.Local),
				time.Date(1989, 8, 31, 23, 59, 59, 0, 0, time.Local))  // 1地域（北見）
		case 2: (time.Date(1989, 9, 26, 0, 0, 0, 0, time.Local),
				time.Date(1989, 6, 4, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 7, 15, 0, 0, 0, 0, time.Local),
				time.Date(1989, 8, 31, 23, 59, 59, 0, 0, time.Local))  // 2地域（岩見沢）
		case 3: (time.Date(1989, 9, 30, 0, 0, 0, 0, time.Local),
				time.Date(1989, 5, 31, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 7, 10, 0, 0, 0, 0, time.Local),
				time.Date(1989, 8, 31, 23, 59, 59, 0, 0, time.Local))  // 3地域（盛岡）
		case 4: (time.Date(1989, 10, 1, 0, 0, 0, 0, time.Local),
				time.Date(1989, 5, 30, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 7, 10, 0, 0, 0, 0, time.Local),
				time.Date(1989, 8, 31, 23, 59, 59, 0, 0, time.Local))  // 4地域（長野）
		case 5: (time.Date(1989, 10, 10, 0, 0, 0, 0, time.Local),
				time.Date(1989, 5, 15, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 7, 6, 0, 0, 0, 0, time.Local),
				time.Date(1989, 8, 31, 23, 59, 59, 0, 0, time.Local))  // 5地域（宇都宮）
		case 6: (time.Date(1989, 11, 4, 0, 0, 0, 0, time.Local),
				time.Date(1989, 4, 21, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 5, 30, 0, 0, 0, 0, time.Local),
				time.Date(1989, 9, 23, 23, 59, 59, 0, 0, time.Local))  // 6地域（岡山）
		case 7: (time.Date(1989, 11, 26, 0, 0, 0, 0, time.Local),
				time.Date(1989, 3, 27, 23, 59, 59, 0, 0, time.Local),
				time.Date(1989, 5, 15, 0, 0, 0, 0, time.Local),
				time.Date(1989, 10, 13, 23, 59, 59, 0, 0, time.Local))  // 7地域（宮崎）
		case 8: (nil, nil,
				time.Date(1989, 3, 25, 0, 0, 0, 0, time.Local),
				time.Date(1989, 12, 14, 23, 59, 59, 0, 0, time.Local))  // 8地域（那覇）
	}

	return (heating_st, heating_en, cooling_st, cooling_en)
}
