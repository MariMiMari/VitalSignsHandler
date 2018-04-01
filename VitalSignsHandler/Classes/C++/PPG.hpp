#include <vector>
#include <math.h>
#include "stats.h"
#include <stdint.h>

#define R_PPG_delay 24

#define PPG_L 120
#define PPG_buf_size 5
#define L2 120
#define PPG_Fs 50 //Hz
#define PPG_filter_order 20
#define PPG_err_lim 100                   //Timeout threshold
#define SpO2_buf_size 3
#define SpO2_cnt_lim 10

typedef struct {	
	bool fault;
	double filtered_RED;
	double filtered_IR;
	double SpO2;
	bool is_new_R;
	double new_R;
	double last_interv;
} PPG_temp;

typedef struct {	
	double RealtimeHR;	// RR по 2 последним пикам(уд/мин)
	double AvgHR;		// RR среднее по выборке(уд/мин)
	double RealtimeRR;	// RR по 2 последним пикам(мс)
	double AvgRR;		// RR среднее по выборке(мс)
	double SpO2;        // SpO2 насыщение кислородом
	double HRV;			//	Изменичвость RR(мс)
	double Filt_RED;	//Последнее число отфильтрованного сигнала красного канала для отрисовки
	double Filt_IR;		//Последнее число отфильтрованного сигнала ИК канала для отрисовки
	int R_occured;
} SPO2_data;

class PPG_class
{
	uint16_t PPG_err_cnt = 0;        //Timeout counter
	bool PPG_R_occured = false;
	double PPG_last_hr[2];

	// HR from PPG
	double RealHR_PPG = 0;
	double AverHR_PPG = 0;
	double SpO2 = 85;

	//	SpO2
	int PPG_SpO2_size;
	vector<double> PPG_SpO2_buf;
	double prev_SpO2_value;
	// Counters of SpO2 timout
	double SpO2_buf[SpO2_buf_size];
	int SpO2_cnt;
	double T2;

	//	 HR
	int PPG_R_size;	vector<double> PPG_R_buf;
	vector<double> PPG_R_del;
	vector<double> PPG_R_time;

	int skip2;
	int L_PPG;
	vector<double> PPG_glob_Red;
	vector<double> PPG_glob_IR;
	bool peak_quit_wind_2;

	int skip = 0;

	// HR from PPG
	double RealHR = 0;
	double AverHR = 0;
	double RealRR = 0;
	double AverRR = 0;
	double HRV = 0;

	//double PPG_signal_RED[PPG_L];
	//double PPG_signal_IR[PPG_L];
	double filt_high_RED[PPG_L];
	double filt_high_IR[PPG_L];
	double Smooth_RED[PPG_L];
	double Smooth_IR[PPG_L];

	// Buffers of input data with fixed sizes
	double PPG_RED_glob[PPG_L];
	double PPG_IR_glob[PPG_L];
	// Buffers for average HR calculation from PPG
	double PPG_interv_buf[PPG_buf_size];
	double PPG_peak_buf[PPG_buf_size];

	PPG_temp My_PPG;	//get_PPG
	SPO2_data Init_PPG_dat(void);
	PPG_temp Init_PPG_Temp(void);
	PPG_temp get_PPG(double *PPG_sign_Red, double *PPG_sign_IR, int PPG_size, double Fs2, int Fconf2, int Skip2);	
public:
	PPG_class();
	SPO2_data glob_res_PPG;	// update_PPG
	SPO2_data update_PPG(double Red, double IR);
};
