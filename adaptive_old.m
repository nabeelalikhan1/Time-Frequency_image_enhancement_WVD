%%%%%%%%%% Main Program for the Adaptive Optimal Kernel %%%%%%%%%%%
% Variable 'sig_in' is the input signal.
% Variable 'ofp' is the output signal.

clear all
addpath('D:\work\tfsa_5-5\windows\win64_bin');
load seizure;
load background;
load S;
% load chirp
% sig_in = chirp;

% load papertest;
% sig_in = papertest;
samp_freq = 512; %in Hz
signal_duration = 0.5;  %in seconds
t = 0:(1/samp_freq):signal_duration-(1/samp_freq);
%sig_in_tmp = cos(2*pi*37*t) + sin(2*pi*63*t) + sin(2*pi*100*t); % 37 + 63 + 100 Hz sinewave
sig_in_tmp = cos(2*pi*120*t.*t+2*pi*90*t) +1* cos(2*pi*120*t.*t+2*pi*60*t) + 0.7*sin(2*pi*40*t)+0.7*sin(2*pi*10*t)+2*exp(-(25*(t-0.35)).^2).*sin(2*pi*72*t);%
%+2*exp(-(25*(t-0.35)).^2).*sin(2*pi*30*t*100);% 37 + 63 + 100 Hz sinewave
sig_in_tmp =2*exp(-(25*(t-0.1)).^2).*sin(2*pi*72*t) +2*exp(-(25*(t-0.35)).^2).*sin(2*pi*72*t);%
t = 1:256;

sig_in_tmp=1*cos(2*pi*(0.48*t-0.0001*t.^2))+cos(6*2*pi.*(cos(2*pi*0.004*t))+2*pi*t*0.2)+cos(6*2*pi.*(-cos(2*pi*0.004*t))+2*pi*t*0.2);%+1*cos(2*pi*(0.05*t+2*0.0000015*t.^2));
%sig_in_tmp=awgn(sig_in_tmp,10,'measured');

%sig_in_tmp = 2*cos(2*pi*-200*t.*t+2*pi*225*t) +2* cos(2*pi*200*t.*t+2*pi*10*t) + 2*sin(2*pi*60*t)+ 2*sin(2*pi*160*t);
%sig_in_tmp=exp(-(25*(t-0.35)).^2).*sin(2*pi*30*t*100);
%sig_in_tmp = cos(2*pi*200*t.*t.*t+2*pi*80*t) +1* cos(2*pi*200*t.*t.*t+2*pi*50*t) + cos(2*pi*200*t.*t.*t+2*pi*20*t);
 %sig_in_tmp=sez_dat(22,:);
% 
%t=t-0.25;
%sig_in_tmp = cos(5*pi*150*t.*t.*t+2*pi*80*t)+cos(5*pi*150*t.*t.*t+2*pi*30*t); 
 load N.mat  % 50 non-seizure segments
% 
% sig_in_tmp = cos(2*pi*100*t.*t+2*pi*100*t) +1* cos(2*pi*100*t.*t+2*pi*130*t) + 0.5*sin(2*pi*5*t)+0.5*sin(2*pi*30*t);
% n=1:256;
% sig_in_tmp = cos(2*pi*0.0004*(n.^2)+2*pi*0.2*n) +1* cos(2*pi*0.0004*(n.^2)+2*pi*0.25*n) + 0.5*sin(2*pi*0.01*(n))+0.5*sin(2*pi*0.06*(n));
fs=samp_freq;

%n = 1:128;

%sig_in_tmp=cos(2*pi*(0.46*n-0.000007*n.^3))+cos(2*pi*(0.4*n-0.000007*n.^3));
%sig_in_tmp=seiz_ds_05_20Hz(200,:);
%sig_in_tmp=sez_dat(200,:);


%sig_in_tmp=back_data(190,:);

%sig_in_tmp = [diff(sig_in_tmp) 0];
sig_in_tmp = hilbert(sig_in_tmp');
% % 
sig_in = [real(sig_in_tmp) imag(sig_in_tmp)];
% 
% 
sig_in_tmp = hilbert(sig_in_tmp');

if ndims(sig_in) ~= 2
    error('I need a 2 x N matrix, real and imaginary parts of the signal.');
end

if size(sig_in, 1) > size(sig_in, 2)
    xr_tmp = sig_in(:,1)'; % IMPORTANT: Pass vectors as rows
    xi_tmp = sig_in(:,2)';  % IMPORTANT: Pass vectors as rows
else
    xr_tmp = sig_in(1, :); % IMPORTANT: Pass vectors as rows
    xi_tmp = sig_in(2, :);  % IMPORTANT: Pass vectors as rows
end

disp('ADAPTIVE OPTIMAL-KERNEL (AOK) TIME-FREQUENCY REPRESENTATION');
disp('   Version 5.0');
disp('  '); % print a blank line
xlen = length(sig_in);

tlag = input('Length of sliding analysis window (power of two, no larger than 256)\n (Number of samples along each dimension of the STAF): ');
disp('  '); % print a blank line
% Check to see if tlag is a power of 2
if (log2(tlag)  - floor(log2(tlag))) ~= 0
    error('Length of sliding window must be a power of 2');
end

% Check to see if tlag is <= 256
if (tlag > 256)
    error('Length of sliding window must be <= 256.');
end

fftlen = input('Number of output frequency samples per time-slice (power of two): ');
disp('  '); % print a blank line

% Check to see if fftlen is a power of 2
if (log2(fftlen) - floor(log2(fftlen))) ~= 0
    error('Number of output frequency samples must be a power of 2');
end

tstep = input('Time increment in samples between time-slice outputs: ');
disp('  '); % print a blank line

if (tstep < 1) || (tstep > xlen)
    error('Time increment must be between 1 and the length of the signal.');
end

% Allocate the output matrix
ofp = zeros(ceil((xlen + tlag + 2)/tstep), fftlen);

vol = input('Normalized volume of optimal kernel (Typically between 1 and 5): ');
disp('  '); % print a blank line

if (vol < 1) || (vol > 5)
    error('Normalized volume of optimal kernel must be between 1 and 5.');
end

% ===========================================================

if (fftlen < (2*tlag))
    fstep = 2*tlag/fftlen;
    fftlen = 2*tlag;
else
    fstep = 1;
end

nits = log2(tstep+2);   % number of gradient steps to take each time

alpha = 0.01;

mu = 0.5;               % gradient descent factor

forget = 0.001;		% set no. samples to 0.5 weight on running AF
nraf = tlag;		% theta size of rectangular AF
nrad = tlag;		% number of radial samples in polar AF
nphi = tlag;		% number of angular samples in polar AF
outdelay = tlag/2;	% delay in effective output time in samples
% nlag-1 < outdelay < nraf to prevent "echo" effect

nlag = tlag + 1;	% one-sided number of AF lags
mfft = ceil(log2(fftlen));
slen = floor( (sqrt(2)) *(nlag-1) + nraf + 3);	% number of delayed samples to maintain
vol = (2.0*vol*nphi*nrad*nrad)/(pi*tlag);   % normalize volume

polafm2 = zeros(nphi, nrad);
rectafr = zeros(nraf, nlag);
rectafi = zeros(nraf, nlag);

xr = zeros(slen,1);
xi = zeros(slen,1);
sigma = ones(nphi,1);

tlen = xlen + nraf + 2;

% ======================================

[rar,rai,rarN,raiN] = rectamake(nlag,nraf,forget);  % make running rect AF parms

plag = plagmake(nrad,nphi,nlag);

[ptheta, maxrad] = pthetamake(nrad,nphi,nraf);   % make running polar AF parms

[rectrotr, rectroti] = rectrotmake(nraf,nlag,outdelay);

[req, pheq] = rectopol(nraf,nlag,nrad,nphi);

outct = 0;

lastsigma = ones(1,nphi);

for ii=0:(tlen-1)

    xr = zeros(1, slen);
    xi = zeros(1, slen);

    if (ii < xlen)
        xr(1:(ii+1)) = fliplr(xr_tmp(1:(ii+1)));
        xi(1:(ii+1)) = fliplr(xi_tmp(1:(ii+1)));
    else
        xr((ii - xlen + 2):(ii + 1)) = fliplr(xr_tmp);
        xi((ii - xlen + 2):(ii + 1)) = fliplr(xi_tmp);
    end

    [rectafr, rectafi] = rectaf(xr,xi,nlag,nraf,rar,rai,rarN,raiN,rectafr,rectafi);

    if ( rem(ii, tstep) == 0 )	% output t-f slice
        outct = outct + 1;

        rectafm2 = rectafr.^2 + rectafi.^2;

        polafm2 = polafint(nrad,nphi,nraf,maxrad,nlag,plag,ptheta,rectafm2);

        %sigma keeps getting updated. need to pass old value into
        %new one

        sigma = sigupdate(nrad,nphi,nits,vol,mu,maxrad,polafm2,lastsigma);
        lastsigma = sigma;

        tfslicer = zeros(1, fftlen);
        tfslicei = zeros(1, fftlen);

        rtemp  = rectafr .* rectrotr + rectafi .* rectroti;
        rtemp1 = rectafi .* rectrotr - rectafr .* rectroti;

        rtemp2 = rectkern(0:(nlag-2),0:(nraf-1),nraf,nphi,req,pheq,sigma);

        % Computer first half of time-frequency domain
        tfslicer(1:(nlag-1)) = sum(rtemp(:, 1:(nlag-1)).*rtemp2);
        tfslicei(1:(nlag-1)) = sum(rtemp1(:, 1:(nlag-1)).*rtemp2);

        % Second half of TF domain is the first half reversed
        tfslicer((fftlen-nlag+3):(fftlen+1)) = tfslicer((nlag-1):-1:1);
        tfslicei((fftlen-nlag+3):(fftlen+1)) = -tfslicei((nlag-1):-1:1);

        % Compute the fft
        % It'd be nice if we could use MATLAB's fft, but I think the
        % custom fft_tfr does some sort of array flipping too
        [tfslicer, tfslicei] = fft_tfr(fftlen,mfft,tfslicer,tfslicei);

        itemp = fftlen/2 + fstep;
        j = 1;
        for i=itemp:fstep:(fftlen-1),
            ofp(outct,j) = tfslicer(i+1);
            j = j + 1;
        end

        for i=0:fstep:(itemp-1),
            ofp(outct,j) = tfslicer(i+1);
            j = j + 1;
        end
    end
end

% Now print the output
f_axis = samp_freq * ((-fftlen/2):fstep:((fftlen/2)-fstep)) / fftlen;  % in Hz
t_axis = 0:tstep:(tlen-1);  % in seconds
contour(t_axis,f_axis,flipud(abs(ofp')));
if length(sig_in_tmp)>200
ofp1=ofp(130:end-129,:);

ofp1=ofp1';
ofp2=ofp1(end/2:end,:);
else

ofp2=ofp(66:end-65,end/2+1:end)';
ofp2=imresize(ofp2,[128 128]);
end


%ofp2=[ofp1(end/2+1:end,:) ofp1(1:end/2,:)];
%  For Gaussian atom
% I1=HTFD_new1(sig_in_tmp,6,15,70);
% I2=HTFD_new1(sig_in_tmp,4,15,70);
% I3=HTFD_new1(sig_in_tmp,3,15,70);
% I4=HTFD_new1(sig_in_tmp,3,20,70);
% I5=HTFD_new1(sig_in_tmp,2,25,70);
% I6=HTFD_new1(sig_in_tmp,1.5,25,70);
% 
% Imin=min(min(I1,I2),min(I3,I4));
% Imin=min(Imin,min(I5,I6));
%mesh(Imin);
%I=HTFD_new(sig_in_tmp,2,5,64);
%I=HTFD_new1(sig_in_tmp,3,35,70*2);
%I1=HTFD_new1(sig_in_tmp,2,25,70*2);


% tfd3 = quadtfd( sig_in_tmp, length(sig_in_tmp)-1, 1, 'specx', 75, 'hamm',256);
% 
% figure;
% tfsapl(sig_in_tmp,tfd3);



[ am I]=wvd1(sig_in_tmp,length(sig_in_tmp));

%%%% 150
C=1;

D=0.065;
%D=0.1;
%%% 173
%D=0.12;

%%% 173
%D=0.12;
g=cskabedbelchourini(length(sig_in_tmp),C,D, D);
am=am.*g;
I = (fft(ifftshift(am,1), [], 1));
I=  ifft(fftshift(I,2), [], 2);
tfd4=real(I);
tfd4(tfd4<0)=0;
figure;%tfsapl(sig_in_tmp,(tfd4));
tfsapl(sig_in_tmp,(tfd4), fs/2,1/fs ,'TimePlot', 'of', 'FreqPlot', 'of','XLabel', 'Frequency(Hz)', 'YLabel', 'Time(s)');


%%%%%%%%%%%%%%%%code fore extended modified B-distribution
wvdz=wvd(sig_in_tmp,length(sig_in_tmp)-1,1,2^nextpow2(length(sig_in_tmp)));
gsig=ifft(fft(wvdz.').');

g=extnd_mbd(0.02,0.2,0.5,length(sig_in_tmp));
smg=gsig.*g;
tfd2=real(fft(ifft(smg.').'));

tfd2(tfd2<0)=0;
figure;%tfsapl(sig_in_tmp,tfd2);
tfsapl(sig_in_tmp,(tfd2), fs/2,1/fs ,'TimePlot', 'of', 'FreqPlot', 'of','XLabel', 'Frequency(Hz)', 'YLabel', 'Time(s)');
% 
% 
 %I=min(I,I1);
% I=min(I,I2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   Directional TFD%%%%%%%%%%%%%%%%%
%I=HTFD_new1(sig_in_tmp,2.5,12,81);
I=DTFD_new(sig_in_tmp,3,12,60);

%tfd_measure(I)
%I(I<0)=0;
figure;%tfsapl(sig_in_tmp,I);
tfsapl(sig_in_tmp,(I), fs/2,1/fs ,'TimePlot', 'of', 'FreqPlot', 'of','XLabel', 'Frequency(Hz)', 'YLabel', 'Time(s)');

figure;%tfsapl(sig_in_tmp,ofp2);
tfsapl(sig_in_tmp,(ofp2), fs/2,1/fs ,'TimePlot', 'of', 'FreqPlot', 'of','XLabel', 'Frequency(Hz)', 'YLabel', 'Time(s)');

% xlabel('Time');
% ylabel('Frequency [Hz]');

tfd=quadtfd(sig_in_tmp,length(sig_in_tmp)-1,1,'specx',55,'hamm');
% fprintf('\nFinished. Output is in variable "ofp"\n');
% figure;contour(t_axis(1:256),f_axis(1:256),flipud(abs(I)));
% figure;contour(t_axis(1:256),f_axis(1:129),flipud(abs(ofp2)));

