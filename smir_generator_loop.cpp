/*
function [ h, H, beta_hat ] = new_smir_generator_with_combined_loop(c, procFs, 
sphLocation, s, L, beta, sphType, sphRadius, mic, N_harm, nsample, K, order, 
refl_coeff_ang_dep, src_ang, src_type)

Inputs:
     c                  speed of sound in m/s
     procFs             processing sampling frequency in Hz
     sphLocation        1 x 3 vector specifying the (x,y,z) coordinates of the 
                        centre of the array in m
     s                  1 x 3 vector specifying the (x,y,z) coordinates of the 
                        source in m
     L                  1 x 3 vector specifying the room dimensions (x,y,z) in m
     beta               1 x 6 vector containing the reflection coefficients or the 
                        "effective" flow resistivity as:
                        [beta_x2 beta_x1 beta_y2 beta_y1 beta_z2 beta_z1] or 
                        beta = Reverberation Time T_60 in s 
     sphType            type of spherical microphone array ['open', 'rigid']
     sphRadius          radius of the spherical microphone array in m
     mic                M x 2 matrix specifying the angles of the microphones
                        (azimuth,inclination) in radians
     N_harm             maximum spherical harmonic order to use in spherical
                        harmonic decomposition
     K                  oversampling factor
     nsample            number of samples of the RIR to calculate 
                        (default=T60*procFs)
     order              reflection order (default=-1, maximum reflection order)
     refl_coeff_ang_dep 0/1; 0 corresponds to real reflection coefficients,
                        1 correspons to angle dependent reflection coefficients
     src_ang            angle of the source in spherical coordinates
     src_type           omnidirectional/subcardioid/cardioid/hypercardioid/bidirectional
 
 Outputs:
     h                  M x nsample matrix containing the calculated RIR(s)
     H                  M x K*nsample/2+1 matrix containing the calculated RTF(s)
     beta_hat           If beta is the reverberation time, the calculated
                        reflection coefficient is returned.
 
 References:
     - D. P. Jarrett, E. A. P. Habets, M. R. P. Thomas, P. A. Naylor,
       "Simulating room impulse responses for spherical microphone arrays,"
       in Proc. IEEE Intl. Conf. on Acoustics, Speech and Signal 
       Processing (ICASSP), May. 2011, pp. 129-132.
     - E. G. Williams, Fourier acoustics: sound radiation and nearfield 
       acoustical holography, 1st ed.	Academic Press, 1999.
     - E. Fisher and B. Rafaely, "The nearfield spherical microphone 
       array," in Proc. IEEE Intl. Conf. on Acoustics, Speech and Signal 
       Processing (ICASSP), Mar. 2008, pp. 5272-5275.
     - J. B. Allen and D. A. Berkley, "Image method for efficiently 
       simulating small-room acoustics", J. Acoust. Soc. Am., vol. 65, 
       no. 4, pp. 943-950, Apr. 1979.
     - Boris Gourevitch and Romain Brette "The impact of early reflections
       on binaural cues", J. Acoust. Soc. Am. Volume 132, Issue 1, 
       pp. 9-27 (2012
     - Takeshi Komatsu, "Improvement of the Delany-Bazley and Miki models
       for fibrous sound-absorbing materials", Acoust. Sci. & Tech. 29, 2 (2008)
 
 This code is based on Emanuel Habets' RIR Generator, available at 
 http://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
 
 Version:          2.0.20130830
 
 History:
    1.0.20101017   Initial version (D. Jarrett)
    1.1.20111212   Performance improvements, added MEX function for most 
                   computationally complex operations, added reflection order (D. Jarrett)
	1.2.20120925   Added truncation of time domain RIRs when oversampling(K > 1) 
                   is used (D. Jarrett)
    2.0.20130830   Main loop in C++ (S. Braun)
                   Added source directivity (S. Braun)
                   Added angle dependent reflection coefficient (S. Braun)
    2.1.20150713   Fixed default RIR length computation

 Copyright (C) 2015 International Audio Laboratories 
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/ 

#define _USE_MATH_DEFINES

#include "iostream"
#include "matrix.h"
#include "mex.h"
#include "math.h"
#include "complex"
#include "vector"
#include "numeric"      // for inner_product
#include "cmath"        // for abs
using namespace std;

complex<double> i = complex<double>(0,1);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

complex<double> refl_factor_komatsu(double phi_img, int N_FFT, double Fs, double beta, int idx);
void normalize(double* arr, int length);
void sphbesselh(int max_nu, double Z, std::vector< complex<double> >& output);
void legendre(int max_n, double X, std::vector<double>& output);
double src_directivity(double* vect1, double* vect2, char src_type);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // Load parameters
    // H = smir_generator_loop_combined(c, procFs, sphLocation, s, L, beta, nsample, order, K, ...
    // shd_k_l_dependent_all_sources, shd_angle_l_dependent_all_sources, mic_pos, sphRadius, ...
    // k, refl_coeff_ang_dep, src_ang, src_type);
    double          c = mxGetScalar(prhs[0]);
    double          fs = mxGetScalar(prhs[1]);
    const double*   rr = mxGetPr(prhs[2]);
    const double*   ss = mxGetPr(prhs[3]);
    const double*   LL = mxGetPr(prhs[4]);
    const double*   beta = mxGetPr(prhs[5]);
    int             nsamples = (int) mxGetScalar(prhs[6]);
    int             order = (int) mxGetScalar(prhs[7]);
    const int       K = (int)mxGetScalar(prhs[8]);
    const double*	shd_k_l_dependent_all_sources_real = mxGetPr(prhs[9]);
    const double*	shd_k_l_dependent_all_sources_imag = mxGetPi(prhs[9]);
    const double*	shd_angle_l_dependent_all_sources = mxGetPr(prhs[10]);
    const double*	mic_pos = mxGetPr(prhs[11]);
    int				M = (int) mxGetM(prhs[11]); // nr of mics
    const double	sphRadius = mxGetScalar(prhs[12]);
    const double*	waveNr = mxGetPr(prhs[13]);
    const int       refl_coeff_ang_dep = (int)mxGetScalar(prhs[14]);
    const double*   src_ang = mxGetPr(prhs[15]); // src_ang in cart coord
    char*           src_type;
    src_type = new char[mxGetN(prhs[16])+1];
    mxGetString(prhs[16], src_type, mxGetN(prhs[16])+1); // int mxGetString(const mxArray *pm, char *str, mwSize strlen);
    
    int k_total = (int) mxGetM(prhs[9]);
    int N_harm = (int) mxGetN(prhs[9])-1;
    double tmp_angle;
    double refl_angles[6];
    int N_FFT = K*nsamples;
    complex<double>  R_p_plus_R_m_beta;
    complex<double> Q[6];
    std::vector<double> legendre_out(N_harm+1);
    std::vector< complex <double> > sphbesselh_out(N_harm+1);
    std::vector< complex <double> > shd_k_l_dependent(k_total*(N_harm+1));
    
    // Create output vector
    plhs[0] = mxCreateDoubleMatrix(M, k_total, mxCOMPLEX);
    double* H_real = mxGetPr(plhs[0]);
    double* H_imag = mxGetPi(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(N_harm+1, M, mxREAL);
    double* shd_angle_l_dependent = mxGetPr(plhs[1]);
    
    // Temporary variables and constants (image-method)
    const double cTs = c/fs;
    double*      r = new double[3];
    double*      s = new double[3];
    double*      L = new double[3];
    double       hu[6];
    double       dist;
    int          fdist;
    int          n1,n2,n3;
    int          q, j, k;
    int          mx, my, mz;
    int          length = 3;
    
    s[0] = ss[0]/cTs; s[1] = ss[1]/cTs; s[2] = ss[2]/cTs;
    L[0] = LL[0]/cTs; L[1] = LL[1]/cTs; L[2] = LL[2]/cTs;
    
    r[0] = rr[0] / cTs;
    r[1] = rr[1] / cTs;
    r[2] = rr[2] / cTs;
    
    n1 = (int) ceil(nsamples/(2*L[0]));
    n2 = (int) ceil(nsamples/(2*L[1]));
    n3 = (int) ceil(nsamples/(2*L[2]));
    
    // Generate room impulse response
    for (mx = -n1 ; mx <= n1 ; mx++)
    {
        hu[0] = 2*mx*L[0];
        
        for (my = -n2 ; my <= n2 ; my++)
        {
            hu[1] = 2*my*L[1];
            
            for (mz = -n3 ; mz <= n3 ; mz++)
            {
                hu[2] = 2*mz*L[2];
                
                for (q = 0 ; q <= 1 ; q++)
                {
                    hu[3] = (1-2*q)*s[0] - r[0] + hu[0];
                    
                    for (j = 0 ; j <= 1 ; j++)
                    {
                        hu[4] = (1-2*j)*s[1] - r[1] + hu[1];
                        
                        for (k = 0 ; k <= 1 ; k++)             
                        {
                            hu[5] = (1-2*k)*s[2] - r[2] + hu[2];
                            
                            dist = sqrt(pow(hu[3], 2) + pow(hu[4], 2) + pow(hu[5], 2)); // = norm(Rp+Rm)/cTs
                            double R_p_plus_R_m[3] = {hu[3]*cTs,hu[4]*cTs,hu[5]*cTs};
                            double R_p_plus_R_m_norm = dist*cTs;
                            
                            if (abs(2*mx-q)+abs(2*my-j)+abs(2*mz-k) <= order || order == -1)
                            {
                                fdist = (int) floor(dist+(sphRadius/cTs));
                                if (fdist < nsamples)
                                {
                                    normalize(R_p_plus_R_m, length);
                                    for (int jj=0; jj<6; jj++){ 
                                        // Calculate angles with the walls
                                        double wall_normal[3] = {0,0,0};
                                        wall_normal[jj/2] = 1;
                                        double dotprod, init = 0.0;
                                        dotprod = std::inner_product(R_p_plus_R_m, R_p_plus_R_m+3, wall_normal, init);
                                        refl_angles[jj] = abs((M_PI/2) - acos(dotprod));
                                    }
                                    
                                    double look_dir_mir[3] = {(pow(-1.0,q))*src_ang[0], (pow(-1.0,j))*src_ang[1], (pow(-1.0,k))*src_ang[2]};     //src_ang in cartesian coord.
                                    double src_rec_vect[3] = {(-1)*R_p_plus_R_m[0],(-1)*R_p_plus_R_m[1],(-1)*R_p_plus_R_m[2]};
                                    
                                    if (sphRadius == 0) {
                                        for (int kk = 0; kk < k_total; kk++) {
                                            complex<double> tmp_H;
                                            if (refl_coeff_ang_dep == 0) {
                                                Q[0] = beta[0]; Q[1] = beta[1]; Q[2] = beta[2];
                                                Q[3] = beta[3]; Q[4] = beta[4]; Q[5] = beta[5];
                                            }
                                            else {
                                                Q[0] = refl_factor_komatsu(refl_angles[0], N_FFT, fs, beta[0], kk); // angle[0] and angle[1] are same but beta[0] and beta[1] may not be same
                                                Q[1] = refl_factor_komatsu(refl_angles[1], N_FFT, fs, beta[1], kk);
                                                Q[2] = refl_factor_komatsu(refl_angles[2], N_FFT, fs, beta[2], kk);
                                                Q[3] = refl_factor_komatsu(refl_angles[3], N_FFT, fs, beta[3], kk);
                                                Q[4] = refl_factor_komatsu(refl_angles[4], N_FFT, fs, beta[4], kk);
                                                Q[5] = refl_factor_komatsu(refl_angles[5], N_FFT, fs, beta[5], kk);
                                            }
                                            R_p_plus_R_m_beta = pow(Q[0],abs(mx-q))*pow(Q[1],abs(mx))*pow(Q[2],abs(my-j))*pow(Q[3],abs(my))*pow(Q[4],abs(mz-k))*pow(Q[5],abs(mz));
                                            R_p_plus_R_m_beta = src_directivity(look_dir_mir, src_rec_vect, src_type[0]) * R_p_plus_R_m_beta;                //attenuation due to directional source
                                            tmp_H =  R_p_plus_R_m_beta * exp(i * waveNr[kk] * R_p_plus_R_m_norm) / R_p_plus_R_m_norm;
                                            for (int ang = 0; ang < M; ang++) {
                                                H_real[ang + M*kk] += tmp_H.real();
                                                H_imag[ang + M*kk] += tmp_H.imag();
                                            }
                                        }
                                    }
                                    else {
                                        for (int ang = 0; ang < M; ang++) {
                                            // Cosine of the angle between the vector R_p+R_m and mic_pos
                                            tmp_angle = R_p_plus_R_m[0] * mic_pos[ang] + R_p_plus_R_m[1] * mic_pos[ang+M] + R_p_plus_R_m[2] * mic_pos[ang+2*M]; // R_p_plus_R_m is normalized
                                            if (tmp_angle < -1)
                                                tmp_angle = -1;
                                            else if (tmp_angle > 1)
                                                tmp_angle = 1;
                                            legendre(N_harm, tmp_angle, legendre_out);
                                            for (int ll = 0; ll <= N_harm; ll++) {  
                                                // Calculating shd_angle_l_dependent; shd_angle_l_dependent(:,l+1) = ((2*l+1) * legendreP(l, 0, tmp_angle)).';
                                                shd_angle_l_dependent[ang + M*ll] = legendre_out[ll] * shd_angle_l_dependent_all_sources[ll];
                                            }
                                        }
                                        
                                        for (int kk = 0; kk < k_total; kk++) {              
                                            // 1i * farfield_mode_strength.' .* repmat(k.',1,N_harm+1) .* besselh(NU1,1,Z1)
                                            sphbesselh(N_harm, waveNr[kk]*R_p_plus_R_m_norm, sphbesselh_out);
                                            for (int ll = 0; ll <= N_harm; ll++) {
                                                shd_k_l_dependent[kk + k_total*ll] = sphbesselh_out[ll] * complex<double>(shd_k_l_dependent_all_sources_real[kk + k_total*ll], shd_k_l_dependent_all_sources_imag[kk + k_total*ll]);
                                            }
                                        }                               
                                        
                                        for (int ang = 0; ang < M; ang++) {
                                            for (int kk = 0; kk < k_total; kk++) {
                                                complex<double> tmp_H;
                                                if (refl_coeff_ang_dep == 0) {
                                                    Q[0] = beta[0]; Q[1] = beta[1]; Q[2] = beta[2];
                                                    Q[3] = beta[3]; Q[4] = beta[4]; Q[5] = beta[5];
                                                }
                                                else {
                                                    Q[0] = refl_factor_komatsu(refl_angles[0], N_FFT, fs, beta[0], kk);
                                                    Q[1] = refl_factor_komatsu(refl_angles[1], N_FFT, fs, beta[1], kk);
                                                    Q[2] = refl_factor_komatsu(refl_angles[2], N_FFT, fs, beta[2], kk);
                                                    Q[3] = refl_factor_komatsu(refl_angles[3], N_FFT, fs, beta[3], kk);
                                                    Q[4] = refl_factor_komatsu(refl_angles[4], N_FFT, fs, beta[4], kk);
                                                    Q[5] = refl_factor_komatsu(refl_angles[5], N_FFT, fs, beta[5], kk);
                                                }
                                                R_p_plus_R_m_beta = pow(Q[0],abs(mx-q))*pow(Q[1],abs(mx))*pow(Q[2],abs(my-j))*pow(Q[3],abs(my))*pow(Q[4],abs(mz-k))*pow(Q[5],abs(mz));
                                                R_p_plus_R_m_beta = src_directivity(look_dir_mir, src_rec_vect, src_type[0]) * R_p_plus_R_m_beta; // Attenuation due to directional source
                                                for (int ll = 0; ll <= N_harm; ll++) {
                                                    tmp_H +=  R_p_plus_R_m_beta * shd_angle_l_dependent[ang + M*ll] * shd_k_l_dependent[kk + k_total*ll];
                                                }
                                                H_real[ang + M*kk] += tmp_H.real();
                                                H_imag[ang + M*kk] += tmp_H.imag();
                                            }
                                        }
                                    }   
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

complex<double> refl_factor_komatsu(double phi_img, int N_FFT, double Fs, double beta, int idx)
{
    double freq;
    complex<double> Z, K, A, B, R;
    if (idx == 0)
        R = complex<double>(1.0,0.0);
    else {
        freq = idx*(Fs/N_FFT);
        Z = 1 + 0.00027*pow((2-log(freq/beta)),6.2) + i*0.0047*pow((2-log(freq/beta)),4.1);         //take conjugate of Z and K
        K = 0.0069*pow((2-log(freq/beta)),4.1) - i*(1 + 0.0004*pow((2-log(freq/beta)),6.2));
        A = sin(phi_img) - pow(Z,-1) * pow(1.0 - pow(K,-2) * pow(cos(phi_img),2),0.5);              //angle(phi) in radians
        B = sin(phi_img) + pow(Z,-1) * pow(1.0 - pow(K,-2) * pow(cos(phi_img),2),0.5);
        R = A/B;
    }
    return R;
}

void normalize(double* arr, int length)
{
    int ii;
    double norm = 0;
    for (ii=0; ii<length; ii++){
        norm += pow(arr[ii],2);
    }
    norm = sqrt(norm);
    for (ii=0; ii<length; ii++){
        arr[ii] = arr[ii]/norm;
    }
}

void sphbesselh(int max_nu, double Z, std::vector< complex<double> >& output)
{
    output[0] = exp(i * Z)/(i * Z);
    output[1] = -i * (-i/Z + 1/(Z*Z)) * exp(i * Z);
    
    for (int nu = 2; nu <= max_nu; nu++) {
        output[nu] = (2*nu-1)/Z * output[nu-1] - output[nu-2];
    }
}

void legendre(int max_n, double X, std::vector<double>& output)
{
    output[0] = 1;
    output[1] = X;
    
    for (int n = 2; n <= max_n; n++) {
        output[n] = (float) (2*n-1)/n * X * output[n-1] - (float) (n-1)/n * output[n-2];
    }
}

double src_directivity(double* vect1, double* vect2, char src_type)
{
    if (src_type=='b' || src_type=='c' || src_type=='s' || src_type=='h')
    {
        double strength, alpha, vartheta, init =0.0;
                
        // Polar Pattern         alpha
        // ---------------------------
        // Bidirectional         0
        // Hypercardioid         0.25
        // Cardioid              0.5
        // Subcardioid           0.75
        // Omnidirectional       1
        
        switch(src_type)
        {
            case 'b':
                alpha = 0;
                break;
            case 'h':
                alpha = 0.25;
                break;
            case 'c':
                alpha = 0.5;
                break;
            case 's':
                alpha = 0.75;
                break;
        };
        
        normalize(vect1, 3);
        normalize(vect2, 3);
        vartheta = std::inner_product(vect1, vect1+3, vect2, init);   //cos(theta)
        strength = alpha + (1-alpha) * vartheta;
        
        return strength;
    }
    else
    {
        return 1;
    }
}