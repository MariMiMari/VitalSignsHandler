#include <vector>
#define M_PI	3.14159265358979323846264338327950288
#define FFT_size 64
#define XYZ_L 128
using namespace std;
double *ComputeLP(int FilterOrder);
double *ComputeHP(int FilterOrder);
double *TrinomialMultiply(int FilterOrder, double *b, double *c);
double *ComputeNumCoeffs(int FilterOrder, double Lcutoff, double Ucutoff, double *DenC);
double *ComputeDenCoeffs(int FilterOrder, double Lcutoff, double Ucutoff);
// �������� � ����� ������� ����� ��������
void add_with_shift(vector<double> &buf, double new_value);
void sort_array(double * arr, int size);
double return_median(double * arr, int size);
//	������� ������������ ����������
double stdev(double *dat, int start_ind, int size);
void add_with_shift(double *buf, double new_value, int size);
//	������� �������������
double min(double *dat, int start_ind, int stop_ind);
//	������� ������������
double max(double *dat, int start_ind, int stop_ind);
//	������� ��������
double mean(vector<double> &dat);
//	������� ��������
double mean(double *dat, int start_ind, int stop_ind);
//	���������� ������ �������. ��������� �������� � ��������� �� ������ ����/2 �� �������������
void Smooth_vector(double *input, double *output, int win_size, int dat_size);
//	���������� ������ �������. ��������� �������� � ��������� �� ������ ����/2 �� �������������
void Smooth_vector(vector<double> &input, vector<double> &output, int win_size);
void FIR_filter(int ord, double *Num, int size, double *x, double *y, double *init);
void IIR_filter(int ord, double *Num, double *Den, int size, double *x, double *y);
//���������� ������ ������ �������� �������������� ��� ��������
double compute_FIR_single(double* input, int filter_order, double *b);
//���������� ������ ������ �������� �������������� ��� ��������
double compute_IIR_single(double* input, double* prev_out, int filter_order, double *b, double *a);
//      ������� �������� ��������� ��������
double non_zero_mean(double *dat, int size);
//	������� ������������������� ��������
double rms(double *dat, int size);
double kurtosis(double *dat, int size);
void fft(double *input, double *out, int size);
void medfilt(double *input, double *output, int arr_size, int win_size);
