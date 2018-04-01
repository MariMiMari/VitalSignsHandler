#include "PPG.hpp"

using namespace std;

PPG_class::PPG_class()
{
    PPG_err_cnt = 0;        //Timeout counter
    PPG_R_occured = false;
    // HR from PPG
    RealHR_PPG = 0;
    AverHR_PPG = 0;
    SpO2 = 85;
    
    //       SpO2
    PPG_SpO2_size = 4;
    PPG_SpO2_buf.resize(PPG_SpO2_size);
    //    prev_SpO2_value;
    // Counters of SpO2 timout
    
    SpO2_cnt = 0;
    T2 = 1.0 / PPG_Fs;
    
    //        HR
    int PPG_R_size = 10;                    //
    PPG_R_buf.resize(PPG_R_size);    // R
    PPG_R_del.resize(PPG_R_size);    // R-R
    PPG_R_time.resize(PPG_R_size);    //
    
    skip2 = 0;    //
    L_PPG = 120;    //
    PPG_glob_Red.resize(L_PPG);    //
    PPG_glob_IR.resize(L_PPG);    //
    peak_quit_wind_2 = true;    //
    
    //    int skip = 0;    //
    
    // HR from PPG
    RealHR = 0;
    AverHR = 0;
    RealRR = 0;
    AverRR = 0;
    HRV = 0;
    
}
//  PPG_data
SPO2_data PPG_class::Init_PPG_dat(void)
{
    SPO2_data temp;
    temp.AvgHR = -1;
    temp.AvgRR = -1;
    temp.HRV = -1;
    temp.RealtimeHR = -1;
    temp.RealtimeRR = -1;
    //temp.size = L_ECG;
    
    temp.Filt_RED = -1;
    temp.Filt_IR = -1;
    
    return temp;
}
//  PPG_temp
PPG_temp PPG_class::Init_PPG_Temp(void)
{
    PPG_temp temp;
    temp.fault = true;
    temp.filtered_RED = 0;
    temp.filtered_IR = 0;
    temp.SpO2 = 0;
    temp.is_new_R = false;
    temp.last_interv = 0;
    temp.new_R = 0;
    return temp;
};
//
PPG_temp PPG_class::get_PPG(double *PPG_sign_Red, double *PPG_sign_IR, int PPG_size, double Fs2, int Fconf2, int Skip2)
{
    //FIR filter coefficiets
    double b_PPG[] = { -0.018460,-0.019277,-0.020068,-0.020831,-0.021564,-0.022265,-0.022931,-0.023561,-0.024152,-0.024703,-0.025212,-0.025678,-0.026098,-0.026472,-0.026799,-0.027077,-0.027305,-0.027484,-0.027612,-0.027689,0.972285,-0.027689,-0.027612,-0.027484,-0.027305,-0.027077,-0.026799,-0.026472,-0.026098,-0.025678,-0.025212,-0.024703,-0.024152,-0.023561,-0.022931,-0.022265,-0.021564,-0.020831,-0.020068,-0.019277,-0.01846 };
    
    PPG_temp Result;    //,
    Result.fault = false;
    
    double RofR;    //
    
    //
    double filt_low_RED = 0;
    double filt_low_IR = 0;
    
    //    double right_shift = 20;     //           R-
    //    int right_offs = 8;    //  ,     R-
    
    //
    if (Fs2 > 25)
    {
        double T = 1.0 / Fs2;    //
        //---------------------------------------------------------------------------------------------------------------
        if (Fconf2 == 2)
        {
            //
            filt_low_RED = mean(PPG_sign_Red, 0, PPG_L);
            filt_low_IR = mean(PPG_sign_IR, 0, PPG_L);
            
            //            double filt_buf[(2 * PPG_filter_order + 1)];
            double x_buf[(2 * PPG_filter_order + 1)];
            
            for (int i = 0; i < (2 * PPG_filter_order + 1); i++)
            {
                // filt_buf[i] = filt_high_RED[PPG_L-(2*PPG_filter_order+1)+i];
                x_buf[i] = PPG_sign_Red[PPG_L - (2 * PPG_filter_order + 1) + i];
            }
            
            add_with_shift(filt_high_RED, compute_FIR_single(x_buf, PPG_filter_order, b_PPG), PPG_L);
            
            for (int i = 0; i < (2 * PPG_filter_order + 1); i++)
            {
                //filt_buf[i] = filt_high_IR[PPG_L-(2*PPG_filter_order+1)+i];
                x_buf[i] = PPG_sign_IR[PPG_L - (2 * PPG_filter_order + 1) + i];
            }
            
            add_with_shift(filt_high_IR, compute_FIR_single(x_buf, PPG_filter_order, b_PPG), PPG_L);
            
            Result.filtered_RED = filt_high_RED[PPG_L - 1];
            Result.filtered_IR = filt_high_IR[PPG_L - 1];
        }
        //----------------------------------------------------------------------------------------------------------------
        
        //       Fconf2
        else
        {
            for (int i = 0; i < PPG_L; i++)
            {
                filt_low_RED = mean(PPG_sign_Red, 0, PPG_L);
                filt_low_IR = mean(PPG_sign_IR, 0, PPG_L);
                for (int i = 0; i < PPG_L; i++)
                {
                    filt_high_RED[i] = PPG_sign_Red[i] - filt_low_RED;
                    filt_high_IR[i] = PPG_sign_IR[i] - filt_low_IR;
                }
            }
        }
        
        //        HR
        double aver_buf[9];
        double aver_buf2[9];
        for (int i = 0; i < 9; i++)
            aver_buf[i] = filt_high_RED[PPG_L - 1 - i];
        Smooth_vector(aver_buf, aver_buf2, 5, 9);
        add_with_shift(Smooth_RED, mean(aver_buf2, 2, 7), PPG_L);
        Result.filtered_RED = Smooth_RED[PPG_L - 1];
        
        for (int i = 0; i < 9; i++)
            aver_buf[i] = filt_high_IR[PPG_L - 1 - i];
        Smooth_vector(aver_buf, aver_buf2, 5, 9);
        add_with_shift(Smooth_IR, mean(aver_buf2, 2, 7), PPG_L);
        Result.filtered_IR = Smooth_IR[PPG_L - 1];
        
        //     5
        //    Smooth_vector(filt_high_RED, Smooth_RED, sm_win_size,PPG_L);
        //    Smooth_vector(filt_high_IR, Smooth_IR, sm_win_size,PPG_L);
        
        //
        double G = max(filt_high_IR, 0, PPG_L) - min(filt_high_IR, 0, PPG_L);
        //   .
        //if (mean(filt_low_IR) > 5000 && G > 25 && G < 1000)
        if (filt_low_IR > 5000 && G > 50 && G < 2000)
        {
            double pks[PPG_L / 2];
            double  locs[PPG_L / 2];
            int N = 0;    //
            //         ( 2     )
            for (int i = PPG_L / 6; i < PPG_L - 3; i++)
            {
                if ((Smooth_IR[i + 1] > Smooth_IR[i]) && (Smooth_IR[i + 1] >= Smooth_IR[i + 2]))
                {
                    pks[N] = Smooth_IR[i + 1];
                    locs[N] = i + 1;
                    N++;
                }
            }
            
            //
            double max_peak = max(filt_high_IR, PPG_L / 6, PPG_L - 2);
            
            double peaks[PPG_L / 4];
            int peaks_ind[PPG_L / 4];
            int K = 0;            //      ,
            double Tres = 0.4;    //           0.3 * .
            double T_delay = 0.3;    //         ()
            double T_200ms = 0.2;
            double temp_delay = 0;
            for (int i = 0; i < N; i++)
            {
                //
                if (pks[i] > Tres*max_peak)
                {
                    //       ,
                    if (K > 0)
                    {    //       ,   T_delay
                        if ((locs[i] - peaks_ind[K - 1])*T > T_200ms)
                        {
                            peaks[K] = pks[i];
                            peaks_ind[K] = locs[i];
                            //max_peak = peaks[K];
                            K = K + 1;
                            //     2
                            //
                            temp_delay = 0.35*(peaks_ind[K - 1] - peaks_ind[K - 2])*T;
                            if (K > 1 && temp_delay > T_200ms)
                            {
                                T_delay = temp_delay;
                            }
                        }
                    }
                    //      ,   ,
                    else
                    {
                        peaks[K] = pks[i];
                        peaks_ind[K] = locs[i];
                        K = K + 1;
                    }
                }
            }
            int Num_of_p = K;
            
            int right_shift = 2;     // offset from right
            int right_offs = 4;    // window width
            
            if (!peak_quit_wind_2)
            {    //        ,
                if (Num_of_p > 1 && ((PPG_L - peaks_ind[Num_of_p - 1]) < (right_shift + right_offs + 1 + Skip2 + right_offs)) && ((PPG_L - peaks_ind[Num_of_p - 1]) > right_shift + right_offs + 1))
                    peak_quit_wind_2 = 1;
            }
            
            //      ,
            if (Num_of_p > 0 && ((PPG_L - peaks_ind[Num_of_p - 1]) < (right_shift + Skip2 + right_offs)) && ((PPG_L - peaks_ind[Num_of_p - 1]) > right_shift))
            {
                if (peak_quit_wind_2)    //    ,
                {
                    Result.is_new_R = 1;   //   -
                    if (Num_of_p > 1)
                        //   -
                        Result.last_interv = peaks_ind[Num_of_p - 1] - peaks_ind[Num_of_p - 2];
                    else
                        Result.last_interv = peaks_ind[Num_of_p - 1];
                    
                    //    -
                    Result.new_R = peaks[Num_of_p - 1];
                }
                else
                {
                    Result.is_new_R = 0;
                    Result.last_interv = 0;
                    Result.new_R = 0;
                }
                peak_quit_wind_2 = 0;    //         ,   ,
            }
            else
            {
                Result.is_new_R = 0;
                Result.last_interv = 0;
                Result.new_R = 0;
            }
            if (Result.is_new_R && K > 1)
            {
                double minRed = Smooth_RED[peaks_ind[K - 2]];
                double minIR = Smooth_IR[peaks_ind[K - 2]];
                for (int i = peaks_ind[K - 2] + 1; i < peaks_ind[K - 1]; i++)
                {
                    if (Smooth_RED[i] < minRed)
                        minRed = Smooth_RED[i];
                    if (Smooth_IR[i] < minIR)
                        minIR = Smooth_IR[i];
                }
                RofR = ((Smooth_RED[peaks_ind[K - 1]] - minRed) / filt_low_RED) / ((Smooth_IR[peaks_ind[K - 1]] - minIR) / filt_low_IR);
                
                //RofR = (rms(filt_high_RED, PPG_L) / filt_low_RED) / (rms(filt_high_IR, PPG_L) / filt_low_IR);
                Result.SpO2 = 108 - 25 * RofR;
            }
            else
                Result.SpO2 = 0;
        }
        //    Kurtosis < 7 ->
        else
        {
            Result.fault = 1;
            Result.is_new_R = 0;
            Result.last_interv = 0;
            Result.new_R = 0;
            Result.SpO2 = 0;
        }
    }
    //        fault   " "
    else
    {
        Result.SpO2 = 0;
        Result.fault = 1;
        Result.is_new_R = 0;
        Result.last_interv = 0;
        Result.new_R = 0;
    }
    //
    return Result;
}
//         1
SPO2_data PPG_class::update_PPG(double Red, double IR)
{
    glob_res_PPG = Init_PPG_dat();
    //    double T = 1 /  PPG_Fs;
    //   /
    double photodiode_ratio = 1.1758;
    //
    //add_with_shift(PPG_glob_Red, Red*photodiode_ratio);
    //add_with_shift(PPG_glob_IR, IR);
    
    add_with_shift(PPG_RED_glob, Red*photodiode_ratio, PPG_L);
    add_with_shift(PPG_IR_glob, IR, PPG_L);
    
    PPG_err_cnt++;
    //add_with_shift(ECG_glob,new_dat,ECG_L);
    // ECG_cnt++;
    My_PPG = get_PPG(PPG_RED_glob, PPG_IR_glob, PPG_L, PPG_Fs, 2, skip);
    //   ,
    if (PPG_err_cnt > PPG_err_lim)
    {
        PPG_last_hr[0] = 1;
        PPG_last_hr[1] = 1;
        PPG_err_cnt = 0;
        for (int i = 0; i < PPG_buf_size; i++)
        {
            PPG_interv_buf[i] = 0;
            PPG_peak_buf[i] = 0;
        }
        for (int i = 0; i < SpO2_buf_size; i++)
        {
            SpO2_buf[i] = 0;
        }
    }
    if (My_PPG.is_new_R != 0)
    {
        PPG_R_occured = true;
        PPG_err_cnt = 0;
        add_with_shift(PPG_last_hr, My_PPG.last_interv, 2);
        if (fabs((60 / (T2*PPG_last_hr[1]) - 60 / (T2*PPG_last_hr[0]))) < 12)
        {
            
            add_with_shift(PPG_interv_buf, My_PPG.last_interv, PPG_buf_size);
            add_with_shift(PPG_peak_buf, My_PPG.new_R, PPG_buf_size);
            if (My_PPG.SpO2 > 100)
                add_with_shift(SpO2_buf, 100, SpO2_buf_size);
            else
                if (My_PPG.SpO2 < 80)
                    add_with_shift(SpO2_buf, 0, SpO2_buf_size);
                else
                    add_with_shift(SpO2_buf, My_PPG.SpO2, SpO2_buf_size); \
        }
    }
    if (non_zero_mean(PPG_interv_buf, PPG_buf_size) == 0)
    {
        AverHR = 0;
        RealHR = 0;
        HRV = 0;
        AverRR = 0;
        RealRR = 0;
    }
    else
    {
        AverHR = 60 / (T2*non_zero_mean(PPG_interv_buf, PPG_buf_size));
        RealHR = 60 / (T2* PPG_interv_buf[PPG_buf_size - 1]);
        AverRR = 1000 * T2*non_zero_mean(PPG_interv_buf, PPG_buf_size);
        RealRR = 1000 * T2*PPG_interv_buf[PPG_buf_size - 1];
        HRV = 1000 * T2*stdev(PPG_interv_buf, 0, PPG_buf_size);
    }
    if (non_zero_mean(SpO2_buf, SpO2_buf_size) == 0)
    {
        glob_res_PPG.SpO2 = 0;
    }
    else
    {
        glob_res_PPG.SpO2 = non_zero_mean(SpO2_buf, SpO2_buf_size);
    }
    glob_res_PPG.Filt_IR = My_PPG.filtered_IR;
    glob_res_PPG.Filt_RED = My_PPG.filtered_RED;
    glob_res_PPG.AvgHR = AverHR;
    glob_res_PPG.RealtimeHR = RealHR;
    glob_res_PPG.AvgRR = AverRR;
    glob_res_PPG.RealtimeRR = RealRR;
    glob_res_PPG.HRV = HRV;
    glob_res_PPG.R_occured = My_PPG.is_new_R;
    return glob_res_PPG;
}
