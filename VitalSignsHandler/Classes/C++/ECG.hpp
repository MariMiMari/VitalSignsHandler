#include <vector>
#include <math.h>
#include "stats.h"
using namespace std;

#define ECG_L 220
#define R_buf_size 10
#define ECG_Fs 100      // Hz
#define filter_conf 2   // Use butter filter
#define ECG_filter_order 4
//-------------------------------Variables for beat classification
#define R_ECG_delay 12    // Delay in samples between is_new_R and beat

#define ECG_err_lim 350      // Timeout threshold
#define classify_size 30

typedef struct {
    double ECG;         // ECG signal
	double RealtimeHR;	// RR по 2 последним пикам(уд/мин)
	double AvgHR;		//* RR среднее по выборке(уд/мин)
	double RealtimeRR;	// RR по 2 последним пикам(мс)
	double AvgRR;		//* RR среднее по выборке(мс)
	double Resp1;		//* оценка частоты дыхани¤(мс)
	double Resp2;
	double HRV;			//* Изменичвость RR(мс)
	double LastFilt;	//Последнее число отфильтрованного сигнала для отрисовки
	double LastSqu;		//Последнее число квадратичного сигнала для отрисовки
	double RespSpline1; //Последнее число сглаженной огибающей амплитуд для отрисовки
	double RespSpline2; //Последнее число сглаженной огибающей интервалов для отрисовки
	int R_occured;
} ECG_data;

typedef struct {
	bool fault;
	double filtered;
	bool is_new_R;
	double new_R;
	double last_interv;
} ECG_temp;

void spline_from_3(double x1, double x2, double x3, double y1, double y2, double y3, double *ab);
//	Возврат интерполированного по сплайну значения
double compute_spline(double *ab, double x1, double x2, double x3, double y1, double y2, double y3, double x);
//	Основная функция работы с ЭКГ
ECG_temp get_ECG(double *ECG_sign, int size, double Fs, int Fconf, int Skip1, double scale);
//	 Добавить во входной массив 1 значение и вычислить параметры ЭКГ
//ECG_data update_ECG(double new_dat, const double Fs, int Fconf, double* beat_stats, int* Pwtt);

class ECG_class
{
	double ECG_glob[ECG_L];		
	double ECG_last_hr[2];
	uint16_t ECG_err_cnt;        //Timeout counter

	double ECG_signal[ECG_L];
	double filtered_ECG[ECG_L];
	//	Возведение сигнала в квадрат
	double squared_ECG[ECG_L];

	// Buffers for average HR calculation from ECG
	double interv_buf[R_buf_size];
	double peak_buf[R_buf_size];

	// HR from ECG
	double RealHR;
	double AverHR;
	double RealRR;
	double AverRR;
	double HRV;
	double T;

	//	Буферы для усредненного вычисления HR и Resp по ЭКГ
	int ECG_R_size;                    //	Размер буферов выходных данных
	vector<double> ECG_R_buf;   //	Буфер для хранения амплитуд R
	vector<double> ECG_R_del;   //	Буфер для хранения интервалов R-R
	vector<double> ECG_R_time;  //	Буфер для хранения абсолютного времени

    //    Счетчик для пропуска вычислений
	int skip;
    
	//Переменная для оценки максимального размаха сигнала ЭКГ
	double Amplitude;

	bool peak_quit_wind;    //	Флаг для отслеживания движения последнего всплеска

	//  Переменные для сглаженного вывода огибающих при оценке дыхания
	int spline_counter;
	int spline_counter_limit;
	vector<double> x_nodes;
	vector<double> y_nodes;

	//	Счетчики и буферы для оценки чаастоты дыхания по сигналу ЭКГ
	int Resp_err_cnt;
	int Resp_err_lim;
	int Resp_buf_size;
	vector<double> Resp1_buf;
	vector<double> Resp2_buf;

	ECG_temp My_ECG; //Экземпляр для возврата из get_ECG

	//	Основная функция работы с ЭКГ
	ECG_temp get_ECG(double *ECG_sign, int size, double Fs, int Fconf, int Skip1, double scale);
	
	//Инициализация экзмемпляра ECG_temp перед началом работы
	ECG_temp Init_Temp(void);
	//Инициализация экзмемпляра ECG_data перед началом работы
	ECG_data Init_ECG_dat(void);
	//	Возврат коэффициентов 3 точечного сплайна
	void spline_from_3(double x1, double x2, double x3, double y1, double y2, double y3, double *ab);
	//	Возврат интерполированного по сплайну значения
	double compute_spline(double *ab, double x1, double x2, double x3, double y1, double y2, double y3, double x);
	
public:
	ECG_data glob_res_ECG;    //Экземпляр для возврата из update_ECG
	ECG_class();	
	//	 Добавить во входной массив 1 значение и вычислить параметры ЭКГ
	ECG_data update_ECG(double new_dat);
};



