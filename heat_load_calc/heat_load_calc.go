package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"time"
)

/*
負荷計算処理の実行

    Args:
        logger
        house_data_path (str): 住宅計算条件JSONファイルへのパス
        output_data_dir (str): 出力フォルダへのパス
        weather_specify_method: 気象データの指定方法
        weather_file_path: 気象データのファイルパス
        region: 地域の区分
*/
func run(
	house_data_path string,
	output_data_dir string,
	weather_specify_method string,
	weather_file_path string,
	region int,
) {
	// interval currently fixed at 15 minutes
	//itv := 15

	// ---- 事前準備 ----

	// 出力ディレクトリの作成
	if _, err := os.Stat(output_data_dir); os.IsNotExist(err) {
		os.Mkdir(output_data_dir, 0755)
	}

	_, err := os.Stat(output_data_dir)
	if os.IsNotExist(err) {
		log.Fatalf("`%s` is not a directory", output_data_dir)
	}

	// 住宅計算条件JSONファイルの読み込み
	log.Printf("住宅計算条件JSONファイルの読み込み開始")
	var rd map[string]interface{}
	if house_data_path[0:4] == "http" {
		resp, err := http.Get(house_data_path)
		if err != nil {
			log.Fatal(err)
		}
		defer resp.Body.Close()
		body, err := ioutil.ReadAll(resp.Body)
		if err != nil {
			log.Fatal(err)
		}
		json.Unmarshal(body, &rd)
	} else {
		file, err := os.Open(house_data_path)
		if err != nil {
			log.Fatal(err)
		}
		defer file.Close()
		bytes, err := ioutil.ReadAll(file)
		if err != nil {
			log.Fatal(err)
		}
		json.Unmarshal(bytes, &rd)
	}

	// 気象データの生成 => weather_for_method_file.csv
	log.Printf("気象データの生成開始")
	w := make_weather(
		weather_specify_method,
		IntervalM15,
		weather_file_path,
		Region(fmt.Sprint(region)),
	)

	log.Printf("スケジュール生成開始")
	scnames := make([]string, 0)
	floor_areas := make([]float64, 0)
	for _, rm := range rd["rooms"].([]interface{}) {
		scnames = append(scnames, rm.(map[string]interface{})["schedule"].(map[string]interface{})["name"].(string))
		floor_areas = append(floor_areas, rm.(map[string]interface{})["floor_area"].(float64))
	}
	scd := get_schedule(
		Auto,
		scnames,
		floor_areas,
	)

	// ---- 計算 ----

	// 計算
	calc(rd, w, scd, IntervalM15, 4, 365, 365, 183)

	// // ---- 計算結果ファイルの保存 ----

	// // 計算結果（瞬時値）
	// result_detail_i_path = path.join(output_data_dir, 'result_detail_i.csv')
	// logger.info('Save calculation results data (detailed version) to `{}`'.format(result_detail_i_path))
	// dd_i.to_csv(result_detail_i_path, encoding='cp932')

	// // 計算結果（平均・積算値）
	// result_detail_a_path = path.join(output_data_dir, 'result_detail_a.csv')
	// logger.info('Save calculation results data (simplified version) to `{}`'.format(result_detail_a_path))
	// dd_a.to_csv(result_detail_a_path, encoding='cp932')
}

func main() {
	var house_data string
	flag.StringVar(&house_data, "input", "", "計算を実行するJSONファイル")

	var output_data_dir string
	flag.StringVar(&output_data_dir, "o", ".", "出力フォルダ")

	var weather string
	flag.StringVar(&weather, "weather", "ees", "気象データの作成方法を指定します。")

	var weather_path string
	flag.StringVar(&weather_path, "weather_path", "", "気象データの絶対パスを指定します。weatherオプションでfileが指定された場合は必ず指定します。")

	var region int
	flag.IntVar(&region, "region", 6, "地域の区分を指定します。気象データの作成方法として建築物省エネ法を指定した場合には必ず指定します。")

	// 引数を受け取る
	flag.Parse()

	if house_data == "" {
		log.Fatal("inputオプションを指定してください。")
	}

	start := time.Now()

	run(
		house_data,
		output_data_dir,
		weather,
		weather_path,
		region,
	)

	elapsedTime := time.Since(start)
	log.Printf("elapsed_time: %v [sec]", elapsedTime)
}
