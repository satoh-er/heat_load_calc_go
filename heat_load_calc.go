package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"time"
)

type Config struct {
	HouseDataPath        string
	OutputDataDir        string
	IsScheduleSaved      bool
	WeatherSpecifyMethod string
	WeatherFilePath      string
	Region               int
	IsWeatherSaved       bool
}

/*
負荷計算処理の実行

    Args:
        logger
        house_data_path (str): 住宅計算条件JSONファイルへのパス
        output_data_dir (str): 出力フォルダへのパス
        is_schedule_saved: スケジュールを出力するか否か
        weather_specify_method: 気象データの指定方法
        weather_file_path: 気象データのファイルパス
        region: 地域の区分
        is_weather_saved: 気象データを出力するか否か
*/
func run(
	house_data_path string,
	output_data_dir string,
	is_schedule_saved bool,
	weather_specify_method string,
	weather_file_path string,
	region int,
	is_weather_saved bool,
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

	// 気象データの保存
	if is_weather_saved {

		weather_path := filepath.Join(output_data_dir, "weather_for_method_file.csv")
		log.Printf("Save weather data to `%s`", weather_path)
		//dd := w.get_weather_as_pandas_data_frame()
		//dd.to_csv(weather_path, encoding='utf-8')
	}

	// スケジュールファイルの保存
	if is_schedule_saved {

		scd.save_schedule(output_data_dir)
	}

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

	var schedule_saved bool
	flag.BoolVar(&schedule_saved, "schedule_saved", false, "スケジュールを出力するか否かを指定します。")

	var weather string
	flag.StringVar(&weather, "weather", "ees", "気象データの作成方法を指定します。")

	var weather_path string
	flag.StringVar(&weather_path, "weather_path", "", "気象データの絶対パスを指定します。weatherオプションでfileが指定された場合は必ず指定します。")

	var region int
	flag.IntVar(&region, "region", 0, "地域の区分を指定します。気象データの作成方法として建築物省エネ法を指定した場合には必ず指定します。")

	var weather_saved bool
	flag.BoolVar(&weather_saved, "weather_saved", false, "気象データを出力するか否かを指定します。")

	var logLevel string
	flag.StringVar(&logLevel, "log", "ERROR", "ログレベルを指定します。 (Default=ERROR)")

	// 引数を受け取る
	flag.Parse()

	// Print flag values
	fmt.Printf("house_data: %s\n", house_data)
	fmt.Printf("output_data_dir: %s\n", output_data_dir)
	fmt.Printf("schedule_saved: %t\n", schedule_saved)
	fmt.Printf("weather: %s\n", weather)
	fmt.Printf("weather_path: %s\n", weather_path)
	fmt.Printf("region: %d\n", region)
	fmt.Printf("weather_saved: %t\n", weather_saved)

	start := time.Now()

	run(
		"example/data_example1.json",
		output_data_dir,
		schedule_saved,
		weather,
		weather_path,
		6,
		weather_saved,
	)

	elapsedTime := time.Since(start)
	log.Printf("elapsed_time: %v [sec]", elapsedTime)
}
