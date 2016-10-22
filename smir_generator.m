function [ h, H, beta_hat ] = smir_generator(c, procFs, sphLocation, s, L, beta, sphType, sphRadius, mic, N_harm, nsample, K, order, varargin)

% function [ h, H, beta_hat ] = smir_generator(c, procFs, sphLocation, s, L, beta, sphType,
%                               sphRadius, mic, N_harm, nsample, K, order, refl_coeff_ang_dep,
%                               HP, src_type, src_ang)
%
% Inputs:
%     c                     speed of sound in m/s
%     procFs                processing sampling frequency in Hz
%     sphLocation           1 x 3 vector specifying the (x,y,z) coordinates of the
%                           centre of the array in m
%     s                     1 x 3 vector specifying the (x,y,z) coordinates of the
%                           source in m
%     L                     1 x 3 vector specifying the room dimensions (x,y,z) in m
%     beta                  1 x 6 vector containing the reflection coefficients or the
%                           "effective" flow resistivity as:
%                           [beta_x2 beta_x1 beta_y2 beta_y1 beta_z2 beta_z1] or
%                           beta = Reverberation Time T_60 in s
%     sphType               type of spherical microphone array ['open', 'rigid']
%     sphRadius             radius of the spherical microphone array in m
%     mic                   M x 2 matrix specifying the angles of the microphones
%                           (azimuth,inclination) in radians
%     N_harm                maximum spherical harmonic order to use in spherical
%                           harmonic decomposition
%     K                     oversampling factor
%     nsample               number of samples of the RIR to calculate
%                           (default=T60*procFs)
%     order                 reflection order (default=-1, maximum reflection order)
%     refl_coeff_ang_dep    0/1; 0 corresponds to real reflection coefficients,
%                           1 correspons to angle dependent reflection coefficients
%     HP                    optional high pass filter (0/1)
%     src_type              omnidirectional/subcardioid/cardioid/hypercardioid/bidirectional
%     src_ang               look angle of the source in spherical coordinates
%
% Outputs:
%     h                     M x nsample matrix containing the calculated RIR(s)
%     H                     M x K*nsample/2+1 matrix containing the calculated RTF(s)
%     beta_hat              If beta is the reverberation time, the calculated
%                           reflection coefficient is returned.
%
% References:
%     - D. P. Jarrett, E. A. P. Habets, M. R. P. Thomas, P. A. Naylor,
%       "Simulating room impulse responses for spherical microphone arrays,"
%       in Proc. IEEE Intl. Conf. on Acoustics, Speech and Signal
%       Processing (ICASSP), May. 2011, pp. 129-132.
%     - E. G. Williams, Fourier acoustics: sound radiation and nearfield
%       acoustical holography, 1st ed.	Academic Press, 1999.
%     - E. Fisher and B. Rafaely, "The nearfield spherical microphone
%       array," in Proc. IEEE Intl. Conf. on Acoustics, Speech and Signal
%       Processing (ICASSP), Mar. 2008, pp. 5272-5275.
%     - J. B. Allen and D. A. Berkley, "Image method for efficiently
%       simulating small-room acoustics", J. Acoust. Soc. Am., vol. 65,
%       no. 4, pp. 943-950, Apr. 1979.
%    -  Boris Gourevitch and Romain Brette "The impact of early reflections
%       on binaural cues", J. Acoust. Soc. Am. Volume 132, Issue 1,
%       pp. 9-27 (2012)
%     - Takeshi Komatsu, "Improvement of the Delany-Bazley and Miki models
%       for fibrous sound-absorbing materials", Acoust. Sci. & Tech. 29, 2 (2008)
%
% This code is based on Emanuel Habets' RIR Generator, available at
% http://home.tiscali.nl/ehabets/rir_generator.html
%
% Version:          2.0.20150713
%  
% History:
%    1.0.20101017   Initial version (D. Jarrett)
%    1.1.20111212   Performance improvements, added MEX function for most 
%                   computationally complex operations, added reflection order (D. Jarrett)
%	 1.2.20120925   Added truncation of time domain RIRs when oversampling(K > 1) 
%                   is used (D. Jarrett)
%    2.0.20130830   Main loop in C++ (S. Braun)
%                   Added source directivity (S. Braun)
%                   Added angle dependent reflection coefficient (S. Braun)
%    2.1.20150713   ixed default RIR length computation
% 
% Copyright (C) 2015 International Audio Laboratories 
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

% error(nargchk(11,13,nargin));

if (nargin < 14)
    refl_coeff_ang_dep = 0;
else
    refl_coeff_ang_dep = varargin{1};
end
if (refl_coeff_ang_dep == 0)
    % Compute wall reflection coefficients from reverberation time
    if (length(beta) == 1)
        V = L(1)*L(2)*L(3);
        S = 2*(L(1)*L(3)+L(2)*L(3)+L(1)*L(2));
        TR = beta;
        % Sabin-Franklin's formula
        alfa = 24*V*log(10)/(c*S*TR);
        if (alfa > 1)
            error('Error: The reflection coefficients cannot be calculated using the room parameters (room dimensions and reverberation time) supplied. Please supply the reflection coefficients or change the reverberation time/room dimensions.');
        end
        beta_hat = sqrt(1-alfa);
        beta = repmat(beta_hat, 1, 6);
    else
        beta_hat = beta;
    end
    
    % Default number of RIR samples to calculate
    if (nargin < 11)
        V = L(1)*L(2)*L(3);
        alpha = ((1-beta(1)^2)+(1-beta(2)^2))*L(2)*L(3) + ...
            ((1-beta(3)^2)+(1-beta(4)^2))*L(1)*L(3) + ...
            ((1-beta(5)^2)+(1-beta(6)^2))*L(1)*L(2);
        TR = 24*log(10.0)*V/(c*alpha);
        if (TR < 0.128)
            TR = 0.128;
        end
        nsample = ceil(TR * procFs);
    end
else
    beta_hat = 0;
end


% Default oversampling factor
if (nargin < 12)
    K = 2;
end
% Default reflection order
if (nargin < 13)
    order = -1;
end

% Default HP
if (nargin < 15)
    HP = 0;
else
    HP = varargin{2};
end

% Default directional src type
if (nargin < 16)
    src_type = 'o';
else
    src_type = varargin{3};
end

% Default look angle towards the receiver
if (nargin < 17)
    [src_ang(1),src_ang(2)] = mycart2sph(sphLocation(1)-s(1),sphLocation(2)-s(2),sphLocation(2)-s(2));
else
    src_ang = varargin{4};
end

% Sanity checks
if (length(sphLocation) ~= 3)
    error('sphLocation must be a 1x3 vector.');
end
if (length(s) ~= 3)
    error('s must be a 1x3 vector.');
end
if (length(L) ~= 3)
    error('L must be a 1x3 vector.');
end
if (~strcmp('open', sphType) && ~strcmp('rigid', sphType))
    error('sphType must be ''open'' or ''rigid''');
end
if (size(mic,2) ~= 2)
    mic = mic.';
    if (size(mic,2) ~= 2)
        error('mic must be an Mx2 matrix.');
    end
end
if (refl_coeff_ang_dep ~= 0 && refl_coeff_ang_dep ~= 1)
    error('refl_coeff_ang_dep must either 0 or 1');
end
if (length(beta) ~= 6 && length(beta) ~= 1)
    error('beta must be a scalar or a 1x6 vector.');
end
if (length(beta) == 1 && refl_coef_ang_dep == 1)
    error('angle dependent reflection coefficients must be a 1x6 vector');
end
if (norm(sphLocation - s) < sphRadius)
    warning('The source cannot be inside the array. No impulse response computed.');
    H = zeros(size(mic, 1), K*nsample/2+1);
    h = 0;
    beta_hat = 0;
    return;
end
if (~all(s <= L) || ~all(s >= 0))
    error('The source must be inside the room.');
end
if (~all(sphLocation + repmat(sphRadius, 1, 3) <= L) || ~all(sphLocation - repmat(sphRadius, 1, 3) >= 0))
    error('The entire array must be inside the room.');
end
if (~strcmp('o', src_type) && ~strcmp('s', src_type) && ~strcmp('c', src_type) && ~strcmp('h', src_type) && ~strcmp('b', src_type))
    error('src_type must be ''o'' or ''s'' or "c" or "h" or "b"');
end
if (size(src_ang) ~= [1,2]) % Check for input parameter 'src_ang'
    src_ang = src_ang.';
    if (size(src_ang) ~= [1,2])
        error('look angle must be a 1x2 matrix.');
    end
end

[src_ang(1), src_ang(2), src_ang(3)] = mysph2cart(src_ang(1),src_ang(2),1);   % Look angle in cartesian coordinates
src_ang = src_ang./norm(src_ang);                                                       % Normalized look angle in cartesian coordinates

N_FFT = K * nsample;                 % Oversampling

[mic_pos(:,1), mic_pos(:,2), mic_pos(:,3)] = mysph2cart(mic(:,1),mic(:,2),1); % Normalised Cartesian coordinates of microphones relative to centre of array (sphLocation)

farfield_mode_strength = zeros(N_harm+1, N_FFT/2+1); % Farfield mode strength - part of the spherical harmonic decomposition (SHD) which is a function of k & l
kk = 1 : N_FFT/2+1;
freq = (kk-1).*procFs/N_FFT;
lambda = c./freq;
k = 2*pi./lambda;
order_l = 0 : N_harm;

% Mode strength calculation
no_overflow = 0;
overflow_warning = 0;
while (no_overflow == 0)
    order_l = 0 : N_harm;
    % Rigid sphere
    if strcmp('rigid', sphType)
        besselh_derivative = repmat(sqrt(pi./(2*k*sphRadius)),N_harm+1,1) .* ( (bsxfun(@besselh,order_l+0.5-1,k'*sphRadius).' - bsxfun(@besselh,order_l+0.5+1,k'*sphRadius).')/2 - bsxfun(@besselh,order_l+0.5,k'*sphRadius).' ./ repmat((2*k*sphRadius),N_harm+1,1) );
        farfield_mode_strength = 1i./(besselh_derivative .* repmat((k*sphRadius).^2,N_harm+1,1)); % Wronskian relation
        
        % If the maximum order N_harm is too high, a bug in MATLAB leads to
        % an overflow in ALL values of besselh (even those corresponding to
        % low orders). This bug been reported to Mathworks and will be
        % fixed in a future release. In the meantime this code will reduce
        % the value of N_harm until no overflow occurs.
        if (any(any(isinf(bsxfun(@besselh,order_l+0.5-1,k'*sphRadius)))))
            % Reduce maximum order N_harm and display warning message
            % (unless it has already been displayed)
            N_harm = N_harm - 1;
            if (overflow_warning == 0)
                warning('The chosen value of N_harm is too high; trying a lower number...');
                overflow_warning = 1;
            end
        else
            no_overflow = 1;
        end
        % Open sphere
    else
        farfield_mode_strength = repmat(sqrt(pi./(2*k*sphRadius)),N_harm+1,1) .* bsxfun(@besselj,order_l+0.5,k'*sphRadius).';
        no_overflow = 1;
    end
end

shd_k_l_dependent_all_sources = 1i * farfield_mode_strength.' .* repmat(k.',1,N_harm+1);
shd_angle_l_dependent_all_sources = (2*order_l+1);

[H] = smir_generator_loop(c, procFs, sphLocation, s, L, beta, nsample, order, K, shd_k_l_dependent_all_sources, shd_angle_l_dependent_all_sources, mic_pos, sphRadius, k, refl_coeff_ang_dep, src_ang, src_type);

H(isnan(H)) = 0;

% To get the RIR the right way around (with MATLAB's Fourier transform
% definition), we need an exponential in the form exp(-ikR)/kR not
% exp(ikR)/kR, so we take the conjugate here. Alternatively Williams
% defines the Fourier transform differently in Fourier Acoustics and
% keeps exp(ikR).
H = conj(H);

h = ifft([H conj(H(:,N_FFT/2:-1:2))], N_FFT, 2, 'symmetric');
% Truncate oversampled RIRs
h = h(:,1:nsample);

if HP==1
    [b,a] = butter(4,50/(procFs/2),'high'); % Cut-off freq is 50 Hz
    h = filter(b,a,h,[],2);
end
