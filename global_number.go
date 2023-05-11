package main

// 空気の比熱, J/kg K
func get_c_a() float64 {
	return 1005.0
}

// 空気の密度, kg/m3
func get_rho_a() float64 {
	return 1.2
}

// 水の蒸発潜熱, J/kg
func get_l_wtr() float64 {
	return 2501000.0
}

// ステファンボルツマン定数
func get_sgm() float64 {
	return 5.67e-8
}

func get_eps() float64 {
	return 0.9
}
