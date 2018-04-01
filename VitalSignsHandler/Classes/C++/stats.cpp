#include "math.h"
#include "stats.h"
#include <algorithm>
#include <complex>

using namespace std;
//    Get min value of input array
double min(double *dat, int start_ind, int stop_ind)
{
    double S = dat[start_ind];
    for (int i = start_ind + 1; i < stop_ind; i++)
    {
        if (dat[i] < S)
            S = dat[i];
    }
    return S;
}

//    Get max value of input array
double max(double *dat, int start_ind, int stop_ind)
{
    double S = dat[start_ind];
    for (int i = start_ind + 1; i < stop_ind; i++)
    {
        if (dat[i] > S)
            S = dat[i];
    }
    return S;
}

//    Get mean value of input vector
double mean(vector<double> &dat)
{
    int size =(int) dat.size();
    double S = 0;
    for (int i = 0; i < size; i++)
    {
        S += dat[i];
    }
    return S / size;
}
//    Get mean value of input array
double mean(double *dat, int start_ind, int stop_ind)
{
    double S = 0;
    for (int i = start_ind; i < stop_ind; i++)
    {
        S += dat[i];
    }
    if (stop_ind - start_ind > 0)
        return S / (stop_ind - start_ind);
    else
        return 0;
}

// Array sort
void sort_array(double* arr, int size)
{
    double temp;
    int i, j;
    for (j = 1; j < size; j++)
    {
        temp = arr[j];
        i = j - 1;
        while (i >= 0 && arr[i] > temp)
        {
            arr[i + 1] = arr[i];
            i = i - 1;
        }
        arr[i + 1] = temp;
    }
}
//Get median value of input arr
double return_median(double* arr, int size)
{
    sort_array(arr, size);
    return arr[size / 2];
}

//median filtering of signal
void medfilt(double *input, double *output, int arr_size, int win_size) {
    //    int i;
    int j;
    double *temp_buf = (double *)malloc(win_size * sizeof(double));
    
    for (j = 0; j < arr_size; j++) {
        if (j <= win_size / 2)
        {
            memcpy(temp_buf, &input[0], sizeof(double)*(j + win_size / 2));
            output[j] = return_median(temp_buf, j + win_size / 2);
        }
        else if (arr_size - j<= win_size / 2)
        {
            memcpy(temp_buf, &input[j - win_size / 2], sizeof(double)*(arr_size - j + win_size / 2));
            output[j] = return_median(temp_buf, (win_size - j) + win_size / 2);
        }
        else
        {
            memcpy(temp_buf, &input[j - win_size / 2], sizeof(double)*(win_size));
            output[j] = return_median(temp_buf, win_size);
        }
    }
    free(temp_buf);
}

//    Get standart deviation
double stdev(double *dat, int start_ind, int size)
{
    double m = mean(dat, start_ind, size);
    double S = 0;
    for (int i = start_ind; i < size; i++)
    {
        S += (dat[i] - m)*(dat[i] - m);//pow((dat[i] - m), 2);
    }
    return sqrt(S / (size - 1 - start_ind));
}

//START------------------------------------------Stuff to compute butter coefficients
double *ComputeLP(int FilterOrder)
{
    double *NumCoeffs;
    int m;
    int i;
    
    NumCoeffs = (double *)calloc(FilterOrder + 1, sizeof(double));
    if (NumCoeffs == NULL) return(NULL);
    
    NumCoeffs[0] = 1;
    NumCoeffs[1] = FilterOrder;
    m = FilterOrder / 2;
    for (i = 2; i <= m; ++i)
    {
        NumCoeffs[i] = (double)(FilterOrder - i + 1)*NumCoeffs[i - 1] / i;
        NumCoeffs[FilterOrder - i] = NumCoeffs[i];
    }
    NumCoeffs[FilterOrder - 1] = FilterOrder;
    NumCoeffs[FilterOrder] = 1;
    
    return NumCoeffs;
}
double *ComputeHP(int FilterOrder)
{
    double *NumCoeffs;
    int i;
    
    NumCoeffs = ComputeLP(FilterOrder);
    if (NumCoeffs == NULL) return(NULL);
    
    for (i = 0; i <= FilterOrder; ++i)
        if (i % 2) NumCoeffs[i] = -NumCoeffs[i];
    
    return NumCoeffs;
}
double *TrinomialMultiply(int FilterOrder, double *b, double *c)
{
    int i, j;
    double *RetVal;
    
    RetVal = (double *)calloc(4 * FilterOrder, sizeof(double));
    if (RetVal == NULL) return(NULL);
    
    RetVal[2] = c[0];
    RetVal[3] = c[1];
    RetVal[0] = b[0];
    RetVal[1] = b[1];
    
    for (i = 1; i < FilterOrder; ++i)
    {
        RetVal[2 * (2 * i + 1)] += c[2 * i] * RetVal[2 * (2 * i - 1)] - c[2 * i + 1] * RetVal[2 * (2 * i - 1) + 1];
        RetVal[2 * (2 * i + 1) + 1] += c[2 * i] * RetVal[2 * (2 * i - 1) + 1] + c[2 * i + 1] * RetVal[2 * (2 * i - 1)];
        
        for (j = 2 * i; j > 1; --j)
        {
            RetVal[2 * j] += b[2 * i] * RetVal[2 * (j - 1)] - b[2 * i + 1] * RetVal[2 * (j - 1) + 1] +
            c[2 * i] * RetVal[2 * (j - 2)] - c[2 * i + 1] * RetVal[2 * (j - 2) + 1];
            RetVal[2 * j + 1] += b[2 * i] * RetVal[2 * (j - 1) + 1] + b[2 * i + 1] * RetVal[2 * (j - 1)] +
            c[2 * i] * RetVal[2 * (j - 2) + 1] + c[2 * i + 1] * RetVal[2 * (j - 2)];
        }
        
        RetVal[2] += b[2 * i] * RetVal[0] - b[2 * i + 1] * RetVal[1] + c[2 * i];
        RetVal[3] += b[2 * i] * RetVal[1] + b[2 * i + 1] * RetVal[0] + c[2 * i + 1];
        RetVal[0] += b[2 * i];
        RetVal[1] += b[2 * i + 1];
    }
    
    return RetVal;
}
double *ComputeNumCoeffs(int FilterOrder, double Lcutoff, double Ucutoff, double *DenC)
{
    double *TCoeffs;
    double *NumCoeffs;
    std::complex<double> *NormalizedKernel;
    double Numbers[11] = { 0,1,2,3,4,5,6,7,8,9,10 };
    int i;
    
    NumCoeffs = (double *)calloc(2 * FilterOrder + 1, sizeof(double));
    if (NumCoeffs == NULL) return(NULL);
    
    NormalizedKernel = (std::complex<double> *)calloc(2 * FilterOrder + 1, sizeof(std::complex<double>));
    if (NormalizedKernel == NULL) return(NULL);
    
    TCoeffs = ComputeHP(FilterOrder);
    if (TCoeffs == NULL) return(NULL);
    
    for (i = 0; i < FilterOrder; ++i)
    {
        NumCoeffs[2 * i] = TCoeffs[i];
        NumCoeffs[2 * i + 1] = 0.0;
    }
    NumCoeffs[2 * FilterOrder] = TCoeffs[FilterOrder];
    double cp[2];
    double Bw, Wn;
    cp[0] = 2 * 2.0*tan(M_PI * Lcutoff / 2.0);
    cp[1] = 2 * 2.0*tan(M_PI * Ucutoff / 2.0);
    
    Bw = cp[1] - cp[0];
    //center frequency
    Wn = sqrt(cp[0] * cp[1]);
    Wn = 2 * atan2(Wn, 4);
    //    double kern;
    const std::complex<double> result = std::complex<double>(-1, 0);
    
    for (int k = 0; k < 11; k++)
    {
        NormalizedKernel[k] = std::exp(-sqrt(result)*Wn*Numbers[k]);
    }
    double b = 0;
    double den = 0;
    for (int d = 0; d < 11; d++)
    {
        b += real(NormalizedKernel[d] * NumCoeffs[d]);
        den += real(NormalizedKernel[d] * DenC[d]);
    }
    for (int c = 0; c < 11; c++)
    {
        NumCoeffs[c] = (NumCoeffs[c] * den) / b;
    }
    
    free(TCoeffs);
    return NumCoeffs;
}
//END------------------------------------------Stuff to compute butter coefficients

double *ComputeDenCoeffs(int FilterOrder, double Lcutoff, double Ucutoff)
{
    int k;            // loop variables
    double theta;     // PI * (Ucutoff - Lcutoff) / 2.0
    double cp;        // cosine of phi
    double st;        // sine of theta
    double ct;        // cosine of theta
    double s2t;       // sine of 2*theta
    double c2t;       // cosine 0f 2*theta
    double *RCoeffs;     // z^-2 coefficients
    double *TCoeffs;     // z^-1 coefficients
    double *DenomCoeffs;     // dk coefficients
    double PoleAngle;      // pole angle
    double SinPoleAngle;     // sine of pole angle
    double CosPoleAngle;     // cosine of pole angle
    double a;         // workspace variables
    
    cp = cos(M_PI * (Ucutoff + Lcutoff) / 2.0);
    theta = M_PI * (Ucutoff - Lcutoff) / 2.0;
    st = sin(theta);
    ct = cos(theta);
    s2t = 2.0*st*ct;        // sine of 2*theta
    c2t = 2.0*ct*ct - 1.0;  // cosine of 2*theta
    
    RCoeffs = (double *)calloc(2 * FilterOrder, sizeof(double));
    TCoeffs = (double *)calloc(2 * FilterOrder, sizeof(double));
    
    for (k = 0; k < FilterOrder; ++k)
    {
        PoleAngle = M_PI * (double)(2 * k + 1) / (double)(2 * FilterOrder);
        SinPoleAngle = sin(PoleAngle);
        CosPoleAngle = cos(PoleAngle);
        a = 1.0 + s2t*SinPoleAngle;
        RCoeffs[2 * k] = c2t / a;
        RCoeffs[2 * k + 1] = s2t*CosPoleAngle / a;
        TCoeffs[2 * k] = -2.0*cp*(ct + st*SinPoleAngle) / a;
        TCoeffs[2 * k + 1] = -2.0*cp*st*CosPoleAngle / a;
    }
    
    DenomCoeffs = TrinomialMultiply(FilterOrder, TCoeffs, RCoeffs);
    free(TCoeffs);
    free(RCoeffs);
    
    DenomCoeffs[1] = DenomCoeffs[0];
    DenomCoeffs[0] = 1.0;
    for (k = 3; k <= 2 * FilterOrder; ++k)
        DenomCoeffs[k] = DenomCoeffs[2 * k - 2];
    
    return DenomCoeffs;
}

// Shift all values of vector to left and add new value to the end
void add_with_shift(vector<double> &buf, double new_value)
{
    rotate(buf.begin(), buf.begin() + 1, buf.end());
    buf.back() = new_value;
}
// Shift all values of array to left and add new value to the end
void add_with_shift(double *buf, double new_value, int size)
{
    for (int i = 0; i < size - 1; i++)
        buf[i] = buf[i + 1];
    buf[size - 1] = new_value;
}

//    Smooth array with average filter of win_size
void Smooth_vector(double *input, double *output, int win_size, int dat_size)
{
    double *temp_buf = (double*)malloc(win_size * sizeof(double));
    for (int i = win_size / 2; i < dat_size - win_size / 2; i++)
    {
        for (int j = 0; j < win_size; j++)
        {
            temp_buf[j] = input[j + i - win_size / 2];
        }
        output[i] = mean(temp_buf, 0, win_size);
    }
    
    for (int i = 0; i < win_size / 2; i++)
    {
        output[i] = input[i];
        output[dat_size - i - 1] = input[dat_size - i - 1];
    }
    free(temp_buf);
}
//    Filter vector with average filter of win_size
void Smooth_vector(vector<double> &input, vector<double> &output, int win_size)
{
    int Length =(int)input.size();
    vector<double> temp_buf(win_size);
    for (int i = win_size / 2; i < Length - win_size / 2; i++)
    {
        for (int j = 0; j < win_size; j++)
        {
            temp_buf[j] = input[j + i - win_size / 2];
        }
        output[i] = mean(temp_buf);
    }
    
    for (int i = 0; i < win_size / 2; i++)
    {
        output[i] = input[i];
        output[Length - i - 1] = input[Length - i - 1];
    }
}
void FIR_filter(int ord, double *Num, int size, double *x, double *y, double *init)
{
    double S;
    for (int i = ord; i < size; i++)
    {
        S = 0;
        for (int j = 0; j < ord + 1; j++)
        {
            S += Num[j] * x[i - j];
        }
        y[i] = S;
    }
    for (int i = 0; i < ord + 1; i++)
        y[i] = init[i];
}

void IIR_filter(int ord, double *Num, double *Den, int size, double *x, double *y)
{
    y[0] = Num[0] * x[0];
    for (int i = 1; i < ord; i++)
    {
        y[i] = 0.0;
        for (int j = 0; j < i; j++)
            y[i] = y[i] + Num[j] * x[i - j];
        for (int j = 1; j < i; j++)
            y[i] = y[i] - Den[j] * y[i - j];
    }
    /* end of initial part */
    for (int i = ord; i < size; i++)
    {
        y[i] = 0;
        for (int j = 0; j < ord; j++)
        {
            y[i] += Num[j] * x[i - j];
        }
        for (int j = 1; j < ord; j++)
        {
            y[i] -= Den[j] * y[i - j];
        }
    }
    
}

// Get new filtered value from input using FIR b coefficients
double compute_FIR_single(double* input, int filter_order, double *b)
{
    double F = 0;
    for (int j = 0; j < (2 * filter_order + 1); j++)
    {
        F = F + b[j] * input[(2 * filter_order) - j];
    }
    return F;
}
// Get new filtered value from input, prev_out using IIR b,a coefficients
double compute_IIR_single(double* input, double* prev_out, int filter_order, double *b, double *a)
{
    double F = 0;
    for (int j = 0; j < (2 * filter_order + 1); j++)
    {
        F = F + b[j] * input[(2 * filter_order) - j];
    }
    for (int j = 1; j < (2 * filter_order + 1); j++)
    {
        F = F - a[j] * prev_out[(2 * filter_order + 1) - j];
    }
    return F;
}

// Get mean value among non-zero items of array
double non_zero_mean(double *dat, int size)
{
    int N = 0;
    double S = 0;
    for (int i = 0; i < size; i++)
    {
        if (dat[i] != 0)
        {
            S += dat[i];
            N++;
        }
    }
    if (N > 0)
        return S / N;
    else
        return 0;
}
//    Get root-mean-square value
double rms(double *dat, int size)
{
    double S = 0;
    for (int i = 0; i < size; i++)
        S += dat[i] * dat[i];
    S /= size;
    return sqrt(S);
}


/**********************************************************************
 FFT - calculates the discrete fourier transform of an array of double
 precision complex numbers using the FFT algorithm.
 
 c = pointer to an array of size 2*N that contains the real and
 imaginary parts of the complex numbers. The even numbered indices contain
 the real parts and the odd numbered indices contain the imaginary parts.
 c[2*k] = real part of kth data point.
 c[2*k+1] = imaginary part of kth data point.
 N = number of data points. The array, c, should contain 2*N elements
 isign = 1 for forward FFT, -1 for inverse FFT.
 */
void FFT(double *c, int N, int isign)
{
    int n, n2, nb, j, k, i0, i1;
    double wr, wi, wrk, wik;
    double d, dr, di, d0r, d0i, d1r, d1i;
    double *cp;
    
    j = 0;
    n2 = N / 2;
    for (k = 0; k < N; ++k)
    {
        if (k < j)
        {
            i0 = k << 1;
            i1 = j << 1;
            dr = c[i0];
            di = c[i0 + 1];
            c[i0] = c[i1];
            c[i0 + 1] = c[i1 + 1];
            c[i1] = dr;
            c[i1 + 1] = di;
        }
        n = N >> 1;
        while ((n >= 2) && (j >= n))
        {
            j -= n;
            n = n >> 1;
        }
        j += n;
    }
    
    for (n = 2; n <= N; n = n << 1)
    {
        wr = cos(2.0 * M_PI / n);
        wi = sin(2.0 * M_PI / n);
        if (isign == 1) wi = -wi;
        cp = c;
        nb = N / n;
        n2 = n >> 1;
        for (j = 0; j < nb; ++j)
        {
            wrk = 1.0;
            wik = 0.0;
            for (k = 0; k < n2; ++k)
            {
                i0 = k << 1;
                i1 = i0 + n;
                d0r = cp[i0];
                d0i = cp[i0 + 1];
                d1r = cp[i1];
                d1i = cp[i1 + 1];
                dr = wrk * d1r - wik * d1i;
                di = wrk * d1i + wik * d1r;
                cp[i0] = d0r + dr;
                cp[i0 + 1] = d0i + di;
                cp[i1] = d0r - dr;
                cp[i1 + 1] = d0i - di;
                d = wrk;
                wrk = wr * wrk - wi * wik;
                wik = wr * wik + wi * d;
            }
            cp += n << 1;
        }
    }
}
//Get amplitude spectrum from input
void fft(double *input, double *out, int size)
{
    double *in;
    in = (double *)malloc(sizeof(double) * size * 2);
    
    for (int i = 0; i < size * 2; i = i + 2)
    {
        in[i] = input[i / 2];
        in[i + 1] = 0;
    }
    
    FFT(in, size, 1);
    for (int i = 0; i < size * 2; i = i + 2)
    {
        out[i / 2] = sqrt(in[i] * in[i] + in[i + 1] * in[i + 1]) / size * 2;
        if (i == 0)
            out[i / 2] /= 2;
    }
}
//Kurtosis of the signal
double kurtosis(double *dat, int size)
{
    double Ku = 0;
    double kUp = 0;
    double kDown = 0;
    float m = mean(dat, 0, size);
    for (int i = 0; i < size; i++)
    {
        kUp += (dat[i] - m)*(dat[i] - m)*(dat[i] - m)*(dat[i] - m);
        kDown += (dat[i] - m)*(dat[i] - m);
    }
    if (kDown != 0)
        Ku = size*kUp / (kDown*kDown);
    return Ku;
}

