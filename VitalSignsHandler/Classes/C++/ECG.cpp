#include "ECG.hpp"
using namespace std;
ECG_class::ECG_class(){
    //Timeout counter
    ECG_err_cnt = 0;
    
    RealHR = 0;
    AverHR = 0;
    RealRR = 0;
    AverRR = 0;
    HRV = 0;
    T = 1.0 / ECG_Fs;
    
    //	Буферы для усредненного вычисления HR и Resp по ЭКГ
    
    //Размер буферов выходных данных
    ECG_R_size = 40;
    //ECG_R_buf.reserve(ECG_R_size);    //	Буфер для хранения амплитуд R
    ECG_R_buf.resize(ECG_R_size);
    //    Буфер для хранения интервалов R-R
    ECG_R_del.resize(ECG_R_size);
    //    Буфер для хранения абсолютного времени
    ECG_R_time.resize(ECG_R_size);
    
    //    Счетчик для пропуска вычислений
    skip = 0;
    
    //Переменная для оценки максимального размаха сигнала ЭКГ
    Amplitude = 0.0000001;
    
    //Флаг для отслеживания движения последнего всплеска
    peak_quit_wind = true;
    
    spline_counter = 0;
    spline_counter_limit = 0;
    x_nodes.resize(3);
    y_nodes.resize(3);
    
    //	Счетчики и буферы для оценки чаастоты дыхания по сигналу ЭКГ
    Resp_err_cnt = 0;
    Resp_err_lim = 1000;
    Resp_buf_size = 5;
    Resp1_buf.resize(Resp_buf_size);
    Resp2_buf.resize(Resp_buf_size);
}

double ECG_class::compute_spline(double *ab, double x1, double x2, double x3, double y1, double y2, double y3, double x)
{
    if (x >= x1 && x <= x2) {
        return ab[0] * pow((x - x1), 3) + ab[1] * (x - x1) + y1;
    }
    else if (x >= x2 && x <= x3) {
        return ab[2] * pow((x - x3), 3) + ab[3] * (x - x3) + y3;
    }
    else if (x > x3)
        return y3;
    else
        return y1;
}

void ECG_class::spline_from_3(double x1, double x2, double x3, double y1, double y2, double y3, double *ab)
{    //a1
    ab[0] = (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) /
    (2 * pow((x1 - x2), 2) * (x1 - x3) * (x2 - x3));
    //b1
    ab[1] = (pow(x1, 2) * (y3 - y2) - x1 * (x2 * (-3 * y1 + y2 + 2 * y3) + 3 * x3 * (y1 - y2)) +
             pow(x2, 2) * (y3 - y1) + x2 * x3 * (y2 - y1) + 2 * pow(x3, 2) * (y1 - y2)) /
    (2 * (x1 - x2) * (x1 - x3) * (x2 - x3));
    //a2
    ab[2] = ab[0] * (x2 - x1) / (x2 - x3);
    //b2
    ab[3] = (2 * pow(x1, 2) * (y2 - y3) + x2 * (x1 * (y3 - y2) + x3 * (2 * y1 + y2 - 3 * y3)) +
             3 * x1 * x3 * (y3 - y2) + pow(x2, 2) * (y3 - y1) + pow(x3, 2) * (y2 - y1)) /
    (2 * (x1 - x2) * (x1 - x3) * (x2 - x3));
}

//Инициализация экзмемпляра ECG_temp перед началом работы
ECG_temp ECG_class::Init_Temp(void) {
    ECG_temp temp;
    temp.fault = true;
    temp.filtered = 0;
    temp.is_new_R = false;
    temp.last_interv = 0;
    temp.new_R = 0;
    return temp;
};

//Инициализация экзмемпляра ECG_data перед началом работы
ECG_data ECG_class::Init_ECG_dat(void) {
    ECG_data temp;
    temp.ECG = -1;
    temp.AvgHR = -1;
    temp.AvgRR = -1;
    temp.HRV = -1;
    temp.RealtimeHR = -1;
    temp.RealtimeRR = -1;
    temp.Resp1 = -1;
    temp.Resp2 = -1;
    //temp.size = L_ECG;
    
    temp.LastFilt = -1;
    temp.LastSqu = -1;
    temp.RespSpline1 = 0;
    temp.RespSpline2 = 0;
    return temp;
}

ECG_temp ECG_class::get_ECG(double *ECG_sign, int size, double Fs, int Fconf, int Skip1, double scale) {
    //Filter coefficients
    double a_ECG[] = { 0, -4.4207, 9.1198, -11.6348, 10.1195, -6.1167, 2.4910, -0.6263, 0.0762 };
    double b_ECG[] = { 0.0186, 0, -0.0743, 0, 0.1114, 0, -0.0743, 0, 0.0186 };
    
    ECG_temp Result;    //Структура, которая вернется из функции
    Result.fault = false;
    // Нормирование сигнала до максимальной амплитуды в 4095
    for (int i = 0; i < ECG_L; i++) {
        ECG_signal[i] = ECG_sign[i] / scale * 4095;
    }
    
    // При малых частотах дискретизации проведение правильных вычислений будет затруднено
    if (Fs > 50) {
        double T = 1.0 / Fs;    //	Период между замерами
        double m = mean(ECG_signal, 0, ECG_L);
        //	При Fconf 1  попытка вырезать из входного сигнала гармонику 50Гц
        //---------------------------------------------------------------Применение полосового фильтра Баттерворта
        if (Fconf == 2) {
            double filt_buf[(2 * ECG_filter_order + 1)];
            double x_buf[(2 * ECG_filter_order + 1)];
            
            for (int i = 0; i < (2 * ECG_filter_order + 1); i++) {
                filt_buf[i] = filtered_ECG[ECG_L - (2 * ECG_filter_order + 1) + i];
                x_buf[i] = ECG_signal[ECG_L - (2 * ECG_filter_order + 1) + i];
            }
            add_with_shift(filtered_ECG,
                           compute_IIR_single(x_buf, filt_buf, ECG_filter_order, b_ECG, a_ECG),
                           ECG_L);
            Result.filtered = filtered_ECG[ECG_L - 1];
        }
        else {//-------------------------------Для других значений Fconf свести фильтрацию к удалению среднего значения
            
            for (int i = 0; i < ECG_L; i++)
                filtered_ECG[i] = ECG_signal[i] - m;
        }
        //-------------------------------------Формирование на выход массива последних значений фильтрованного сигнала
        double aver_buf[9];
        double aver_buf2[9];
        for (int i = 0; i < 9; i++)
            aver_buf[i] = filtered_ECG[ECG_L - 1 - i] * filtered_ECG[ECG_L - 1 - i];
        Smooth_vector(aver_buf, aver_buf2, 5, 9);
        add_with_shift(squared_ECG, mean(aver_buf2, 2, 7), ECG_L);
        //Result.filtered = squared_ECG[ECG_L - 1] / -100;
        
        //-----------Вычисления коэффициента эксцесса (или куртозиса - косвенная оценка наличия узких пиков в сигнале)
        double Ku = kurtosis(filtered_ECG, ECG_L);
        
        //Оценка разброса сигнала относительно его среднего значения
        //double G = 2 * stdev(ECG_signal, ECG_L / 5, ECG_L) / abs(mean(ECG_signal, ECG_L / 5, ECG_L));
        double MAX = max(ECG_signal, ECG_L / 4, ECG_L);
        double MIN = min(ECG_signal, ECG_L / 4, ECG_L);
        double G = stdev(ECG_signal, ECG_L / 4, ECG_L) / (MAX - MIN);
        // G>1 -> сигнал слишком зашумлен
        //	Kurtosis > 4  -> после фильтрации остались узкие высокие всплески, с которыми можно дальше работать
        if (Ku > 4 && G < 0.4) {
            static double pks[ECG_L / 2];
            static double locs[ECG_L / 2];
            int pks_N = 0;    //	Счетчик для всех локальных максимумов
            //	Поиск максимумов	(Игнорирование 2 значений в конце и начале)
            for (int i = ECG_L / 5; i < ECG_L - 3; i++) {
                if ((squared_ECG[i + 1] > squared_ECG[i]) &&
                    (squared_ECG[i + 1] >= squared_ECG[i + 2])) {
                    pks[pks_N] = squared_ECG[i + 1];
                    locs[pks_N] = i + 1;
                    pks_N++;
                }
            }
            //   Выбор для первых измерений в качестве компонента порога максимального пика
            double max_peak = max((double *)(squared_ECG), ECG_L / 5, ECG_L - 3);
            
            static double peaks[ECG_L / 6];
            static double peaks_ind[ECG_L / 6];
            int peaks_N = 0;            //	Счетчик для пиков, прошедших через порог
            double Tres = 0.30;    //	Для поиска следующего всплеска на уровне выше 0.3 * Ампл. текущего
            double T_delay = 0.3;    //Минимальная задержка по времени между соседними пиками при первом измерении(с)
            double T_200ms = 0.2;
            double temp_delay = 0;
            for (int i = 0; i < pks_N; i++) {
                // Проверка превышение пиком порога
                if (pks[i] > Tres * max_peak) {
                    //   Если один пик уже есть,
                    if (peaks_N > 0) {    //	искать следующий не раньше, чем черезез T_delay
                        if ((locs[i] - peaks_ind[peaks_N - 1]) * T > T_200ms) {
                            peaks[peaks_N] = pks[i];
                            peaks_ind[peaks_N] = locs[i];
                            max_peak = peaks[peaks_N];
                            peaks_N++;
                            //   После обнаружении 2го пика начать подстройку
                            //   временной задержки
                            temp_delay =
                            0.35 * (peaks_ind[peaks_N - 1] - peaks_ind[peaks_N - 2]) * T;
                            if (peaks_N > 1 && temp_delay > T_200ms) {
                                T_delay = temp_delay;
                            }
                        }
                    }
                    //   Если самый первый всплеск, прошедший через порог, просто записать его
                    else {
                        peaks[peaks_N] = pks[i];
                        peaks_ind[peaks_N] = locs[i];
                        peaks_N++;
                    }
                }
            }
            int Num_of_p = peaks_N;    //Общее число найденных R-всплесков
            int right_shift = 3;
            if (!peak_quit_wind) {    //	Предыдущий всплеск покинул первое окно, поднять флаг
                if (Num_of_p > 0 && ((ECG_L - peaks_ind[peaks_N - 1]) < (right_shift + 7 + 3)) &&
                    ((ECG_L - peaks_ind[peaks_N - 1]) > right_shift + 7))
                    peak_quit_wind = 1;
            }
            
            //	Поиск нового пика, который лежит в окне
            if (Num_of_p > 0 && ((ECG_L - peaks_ind[peaks_N - 1]) < (right_shift + 3)) &&
                ((ECG_L - peaks_ind[peaks_N - 1]) > right_shift)) {
                if (peak_quit_wind)    //	Проверка, что это не предыдущий всплеск
                {
                    Result.is_new_R = 1;   // Первый элемент - обнаружен новый пик
                    if (Num_of_p > 1)
                        // Второй элемент - расстояние до предпоследнего пика
                        Result.last_interv = peaks_ind[peaks_N - 1] - peaks_ind[peaks_N - 2];
                    else
                        Result.last_interv = peaks_ind[peaks_N - 1];
                    
                    //  Третий элемент - значение последнего пика
                    Result.new_R = peaks[peaks_N - 1];
                }
                else {
                    Result.is_new_R = 0;
                    Result.last_interv = 0;
                    Result.new_R = 0;
                }
                
                peak_quit_wind = 0;    //	при обнаружении нового пика опустить флаг, до тех пор, пока он не покинет первое окно
            }
            else {
                Result.is_new_R = 0;
                Result.last_interv = 0;
                Result.new_R = 0;
            }
            
        }
        //	Ku < 4 && G < 1-> скорее всего на входной выборке нет ЭКГ всплесков
        else {
            Result.fault = true;
            Result.is_new_R = false;
            Result.last_interv = 0;
            Result.new_R = 0;
        }
    }
    // Если частота дискретизации слишком мала поднять флаг fault и вернуть "пустую структуру"
    else {
        Result.filtered = 0;
        Result.fault = true;
        Result.is_new_R = false;
        Result.last_interv = 0;
        Result.new_R = 0;
    }
    //	Возврат сформированной структуры
    return Result;
}

ECG_data ECG_class::update_ECG(double new_dat) {
    
    add_with_shift(ECG_glob, new_dat, ECG_L);
    glob_res_ECG = Init_ECG_dat();
    glob_res_ECG.ECG = new_dat;
    //Remember max
    if (max(ECG_glob, 0, ECG_L) > Amplitude && max(ECG_glob, 0, ECG_L) < 65536 + 1)
        Amplitude = max(ECG_glob, 0, ECG_L);
    if ((-1 * min(ECG_glob, 0, ECG_L)) > Amplitude && (-1 * min(ECG_glob, 0, ECG_L)) < 65536 + 1)
        Amplitude = -min(ECG_glob, 0, ECG_L);
    
    ECG_err_cnt++;
    My_ECG = get_ECG(ECG_glob, ECG_L, ECG_Fs, filter_conf, skip, Amplitude);
    
    // Сигнал надолго пропал, обнуление буферов выходных данных
    if (ECG_err_cnt > ECG_err_lim) {
        ECG_err_cnt = 0;
        ECG_last_hr[0] = 1;
        ECG_last_hr[1] = 1;
        for (int i = 0; i < R_buf_size; i++) {
            interv_buf[i] = 0;
            peak_buf[i] = 0;
        }
        for (int i = 0; i < ECG_R_size; i++) {
            ECG_R_buf[i] = 0;
            ECG_R_time[i] = 0;
            ECG_R_del[i] = 0;
        }
    }
    
    if (My_ECG.is_new_R != 0) {
        ECG_err_cnt = 0;
        add_with_shift(ECG_last_hr, My_ECG.last_interv, 2);
        if (abs((int)(60 / (T * ECG_last_hr[1]) - 60 / (T * ECG_last_hr[0]))) < 16) {
            add_with_shift(interv_buf, My_ECG.last_interv, R_buf_size);
            add_with_shift(peak_buf, My_ECG.new_R, R_buf_size);
            
            add_with_shift(ECG_R_buf, My_ECG.new_R);
            add_with_shift(ECG_R_del, My_ECG.last_interv);
            //My_ECG.is_new_R = false;
            ECG_err_cnt = 0;
            spline_counter = 0;
        }
    }
    if (non_zero_mean(interv_buf, R_buf_size) == 0) {
        AverHR = 0;
        RealHR = 0;
        HRV = 0;
        AverRR = 0;
        RealRR = 0;
    }
    else {
        AverHR = 60 / (T * non_zero_mean(interv_buf, R_buf_size));
        RealHR = 60 / (T * interv_buf[R_buf_size - 1]);
        AverRR = 1000 * T * non_zero_mean(interv_buf, R_buf_size);
        RealRR = 1000 * T * interv_buf[R_buf_size - 1];
        HRV = 1000 * T * stdev(interv_buf, 0, R_buf_size);
    }
    
    glob_res_ECG.AvgHR = AverHR;
    glob_res_ECG.AvgRR = AverRR;
    glob_res_ECG.HRV = HRV;
    glob_res_ECG.LastFilt = My_ECG.filtered;
    glob_res_ECG.LastSqu = My_ECG.filtered;
    glob_res_ECG.RealtimeHR = RealHR;
    glob_res_ECG.RealtimeRR = RealRR;
    glob_res_ECG.R_occured = My_ECG.is_new_R;
    
    //    Формирование абсолютного времени
    ECG_R_time[0] = ECG_R_del[0];
    for (int i = 1; i < ECG_R_size; i++) {
        ECG_R_time[i] = ECG_R_time[i - 1] + ECG_R_del[i];
    }
    
    //------------------------------------Двукратное сглаживание пиков и интервалов для оценки RespR
    vector<double> Sm_Resp1(ECG_R_size);
    vector<double> Sm_Resp2(ECG_R_size);
    //vector<double> temp1(ECG_R_size);
    //vector<double> temp2(ECG_R_size);
    Smooth_vector(ECG_R_buf, Sm_Resp1, 3);
    Smooth_vector(ECG_R_del, Sm_Resp2, 3);
    /*for (int i = 0; i < ECG_R_size; i++)
     Sm_Resp2[i] = ECG_R_del[i];*/
    //------------------------------------Двукратное сглаживание пиков и интервалов для оценки RespR
    //Smooth_vector(ECG_R_buf, Sm_Resp1, 3);
    //Smooth_vector(ECG_R_del, Sm_Resp2, 3);
    //Smooth_vector(temp1, Sm_Resp1, 3);
    //Smooth_vector(temp2, Sm_Resp2, 3);
    //---------------------------------------------------------------------------------------------
    
    //----------Получение последнего значения сплайна для отрисовки дыхания по сглаженным амплитудам
    double y;
    x_nodes[0] = 0;
    x_nodes[1] = ECG_R_del[ECG_R_size - 3];
    x_nodes[2] = ECG_R_del[ECG_R_size - 2] + ECG_R_del[ECG_R_size - 3];
    
    y_nodes[0] = Sm_Resp1[ECG_R_size - 4];
    y_nodes[1] = Sm_Resp1[ECG_R_size - 3];
    y_nodes[2] = Sm_Resp1[ECG_R_size - 2];
    
    if (ECG_R_del[ECG_R_size - 6] != 0) {
        double AB[4];
        spline_from_3(x_nodes[0], x_nodes[1], x_nodes[2], y_nodes[0], y_nodes[1], y_nodes[2], AB);
        
        spline_counter_limit = ECG_R_del[ECG_R_size - 3];
        if (spline_counter < spline_counter_limit + 1) {
            y = compute_spline(AB, x_nodes[0], x_nodes[1], x_nodes[2], y_nodes[0], y_nodes[1],
                               y_nodes[2], spline_counter);
        }
        else
            y = y_nodes[1];
    }
    else {
        y = 0;
    }
    glob_res_ECG.RespSpline1 = y;
    //----------------------------------------------------------------------------------------------
    
    //----------Получение последнего значения сплайна для отрисовки дыхания по сглаженным интервалам
    x_nodes[0] = 0;
    x_nodes[1] = ECG_R_del[ECG_R_size - 3];
    x_nodes[2] = ECG_R_del[ECG_R_size - 2] + ECG_R_del[ECG_R_size - 3];
    
    y_nodes[0] = Sm_Resp2[ECG_R_size - 4];
    y_nodes[1] = Sm_Resp2[ECG_R_size - 3];
    y_nodes[2] = Sm_Resp2[ECG_R_size - 2];
    
    if (ECG_R_del[ECG_R_size - 6] != 0) {
        double AB[4];
        spline_from_3(x_nodes[0], x_nodes[1], x_nodes[2], y_nodes[0], y_nodes[1], y_nodes[2], AB);
        
        spline_counter_limit = ECG_R_del[ECG_R_size - 3];
        if (spline_counter < spline_counter_limit + 1) {
            y = compute_spline(AB, x_nodes[0], x_nodes[1], x_nodes[2], y_nodes[0], y_nodes[1], y_nodes[2], spline_counter);
        }
        else
            y = y_nodes[1];
    }
    else {
        y = 0;
    }
    glob_res_ECG.RespSpline2 = y;// / (T_19*Fs);
    
    spline_counter++;
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------Поиск максимумов огибающей кривой и их координат
    vector<double> pks1;
    vector<double> pks2;
    vector<double> locs1;
    vector<double> locs2;
    for (int i = 0; i < ECG_R_size - 2; i++) {
        if (Sm_Resp1[i] < Sm_Resp1[i + 1] && Sm_Resp1[i + 2] <= Sm_Resp1[i + 1]) {
            pks1.push_back(Sm_Resp1[i + 1]);
            locs1.push_back(ECG_R_time[i + 1]);
        }
        if (Sm_Resp2[i] < Sm_Resp2[i + 1] && Sm_Resp2[i + 2] <= Sm_Resp2[i + 1]) {
            pks2.push_back(Sm_Resp2[i + 1]);
            locs2.push_back(ECG_R_time[i + 1]);
        }
    }
    //----------------------------------------------------------------------------------------------
    
    //------------------------------------------------------------Поиск интервалов между максимумами
    if (pks1.size() > 1 && pks2.size() > 1) {
        vector<double> interv1;
        vector<double> interv2;
        for (int i = 1; i < locs1.size(); i++)
            interv1.push_back(locs1[i] - locs1[i - 1]);
        for (int i = 1; i < locs2.size(); i++)
            interv2.push_back(locs2[i] - locs2[i - 1]);
        //------------------------------------------------------------------------------------------
        
        //-----------------------------------------------------Вычисление частоты во вдохах в минуту
        for (int i = 0; i < Resp_buf_size - 1; i++) {
            Resp1_buf[i] = Resp1_buf[i + 1];
            Resp2_buf[i] = Resp2_buf[i + 1];
        }
        Resp1_buf[Resp_buf_size - 1] = 60 / (T * mean(interv1));
        Resp2_buf[Resp_buf_size - 1] = 60 / (T * mean(interv2));
        glob_res_ECG.Resp1 = mean(Resp1_buf);
        glob_res_ECG.Resp2 = mean(Resp2_buf);
        //--------------Ограниечение не чаще раз в 2 секунды
        if (glob_res_ECG.Resp1 > 30)
            glob_res_ECG.Resp1 = 30;
        if (glob_res_ECG.Resp2 > 30)
            glob_res_ECG.Resp2 = 30;
        //------------------------------------------------------------------------------------------
        Resp_err_cnt = 0;
    }
    else {
        Resp_err_cnt++;
        if (Resp_err_cnt > Resp_err_lim)
            for (int i = 0; i < Resp_buf_size - 1; i++) {
                Resp1_buf[i] = 0;
                Resp2_buf[i] = 0;
            }
        glob_res_ECG.Resp1 = 0;
        glob_res_ECG.Resp2 = 0;
    }
    return glob_res_ECG;
}
