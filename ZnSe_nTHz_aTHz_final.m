
% Here I considered two options of calculation of nTHz. 
%1st option: first I calculated aTHz  from equation 3 in book:X.-C. Zhang, J. Xu, Introduction
%to THz Wave Photonics, after with the help of The Kramers-Kronig relation for the real part n(ω), this is refractive index in THz range. 
%2nd option: aTHz was calculated by the same way in 1st method. refractive
%index in THz was calculated by equation 3 of the same book
clear; clc;
cry = 7; % 0 - LN 3 - Gap 4 - GaAs  7 - ZnSe  2 - ZnTe 
files = {'ref 9e3 avg_0', '9e3_0'}; 
ETHz = length(files);

dataCells = cell(ETHz, 1);

for i = 1:ETHz
    dataCells{i} = readmatrix(files{i});
end


Etref = dataCells{1};   % the reference signal in time domain
Etsam = dataCells{2};   % the sample signal in time domain

Eref = fft(Etref);
Esam = fft(Etsam);
phi_ref = angle(Eref);
phi_sam = angle(Esam);
delta_phi = - unwrap(phi_sam - phi_ref); %phase difference
% Esam = abs(fft(Etsam)); %the magnitude spectrum for spectral analysis.
% But for the calculation of nTHz and aTHz I used simple fft(Esam) fft(Etsam)and
% fft(Etref).
%sam = abs(fft(Etsam).^2);
N = length(Etref);
c0 = 3e8;
lambda0 = 10.6e-6;
omega0 = 2*pi*c0/lambda0;
omegaMAX = 5*2*omega0;
domega = omegaMAX/N;
omega = (0:N-1)*domega;
nu=omega/2/pi;
n0 = 1;
d = 1.0e-3;

 
M = floor(N/2); 
T = 300;
dt = T/(N-1);
freq = [0:N-1]'/T * 1e-12;
time=Etref(:,1);

%  Calculate Phase Shift
%Phase_shift = omega * Delta_t;  % Phase shift as a function of frequency

%n = 1 + (c * Phase_shift) / (omega * z_vegso);
%Transmission =  Esam(:, 2) ./ Eref(:, 2);
%n = 1.2321;
%aTHz = (2 ./ z_vegso) .* log((abs(Eref(:, 2)) ./ abs(Esam(:, 2))) .* ((1 + n)^2 / 4 * n));

% max value of Etref
%[max_val_1, time1] = max(Etref(:, 2));  % Find the maximum value in the 2nd column of reference signal in time domain
%max_Etref = Etref(time1, 1);  % corresponding value from the 1st column which is time
%%max value of Etsam
% max value of Etref
%[max_val_2, time2] = max(Etsam(:, 2));  % Find the maximum value in the 2nd column of reference signal in time domain
%max_Etsam = Etsam(time2, 1);  % corresponding value from the 1st column which is time

% Tine delay
%delta_t = max_Etsam - max_Etref;  % in ps
%aTHz = (1 ./ z_vegso) * log(abs(Eref(:, 2)) ./ abs(Esam(:, 2)));    % from X.-C. Zhang, J. Xu, Introduction to THz Wave Photonics,



H1=Esam./Eref;
H = Esam(:, 2) ./ Eref(:, 2);
n = 1 + delta_phi(:,2)' * c0 ./ omega * d;
k1 = log(4*n'./(abs(H).*(n'+1).^2))*c0./freq/2/pi/d; %k1- absorption index
alpha = 2*2*pi*freq.*k1/c0; %absorption coefficient
n_comp = n-i*k1;

subplot(2,2,1);
plot(time,Etref,time,Etsam);
subplot(2,2,2);
plot(freq,abs(Eref),freq,abs(Esam));

n_sg = sgolayfilt(n,2,5);
alpha_sg = sgolayfilt(alpha,2,5);

subplot(2,2,3);
plot(freq,n,freq,n_sg);

axis([0.2 3 3 7]);

subplot(2,2,4);
plot(freq,alpha,freq,alpha_sg);
axis([0.2 3 0 20]);

%% Export data for aTHz
exportDataToFile('aTHz.txt', 'THz Frequency, THz\taTHz\n', freq,alpha);

% Export data for nTHz
exportDataToFile('nTHz.txt', 'THz Frequency, THz\tnTHz\n', freq, n);
exportDataToFile('nTHz_sg.txt', 'THz Frequency, THz\tnTHz\n', freq, n_sg);
exportDataToFile('aTHz_sg.txt', 'THz Frequency, THz\taTHz\n', freq, alpha_sg);
exportDataToFile('Spectrum.txt', 'THz Frequency, THz\Spectrum\n', freq, Esam);
% Define the common export function
function exportDataToFile(filename, headers, col1, col2)
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Check if file opened successfully
    if fileID == -1
        error('Could not open file %s for writing.', filename);
    end

    % Write the headers
    fprintf(fileID, '%s\n', headers);

    % Write data side-by-side
    for i = 1:length(col1)
        fprintf(fileID, '%f\t%f\n', col1(i), col2(i));
    end

    % Close the file
    fclose(fileID);
    disp(['Data exported to ', filename]);
end

% Define data and headers
% nuu = [0.1, 0.2, 0.3];      % Example frequency data (replace with actual data)
% abso = [10, 20, 30];        % Example absorption data for aTHz
% nTHZ = [1.5, 1.6, 1.7];     % Example refractive index data for nTHz




