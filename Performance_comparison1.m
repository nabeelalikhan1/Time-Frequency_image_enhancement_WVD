clear all;
close all;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 1:256;
%t = -127:128;
%load seizure;

tmp_sig=4*exp(-0.001*(t-192).^2).*(cos((2*pi*0.07*t)))+2*cos(2*pi*(0.05*t+0.0000015*t.^3))+2*cos(2*pi*(0.075*t+0.0000015*t.^3))+4*exp(-0.001*(t-64).^2).*(cos((2*pi*0.45*t)))+4*exp(-0.001*(t-192).^2).*(cos((2*pi*0.45*t)));%+1*cos(2*pi*(0.45*t));
%tmp_sig =sez_dat(63,:);
%tmp_sig=2*cos(2*pi*(0.05*t+0.0005*t.^2))+2*cos(2*pi*(0.25*t-0.0005*t.^2))+2*exp(-0.001*(t-64).^2).*(cos((2*pi*0.45*t)))+2*exp(-0.001*(t-192).^2).*(cos((2*pi*0.45*t)));%+1*cos(2*pi*(0.45*t));

%tmp_sig=2*cos(2*pi*(0.05*t+0.00075*t.^2))+2*cos(2*pi*(0.45*t-0.00075*t.^2));

%tmp_sig=2*cos(2*pi*(0.05*t+0.0000015*t.^3))+2*cos(2*pi*(0.15*t+0.0000015*t.^3))+1*cos(2*pi*(0.45*t));

%tmp_sig=awgn(tmp_sig,7,'measured');
% y=rand(1,length(tmp_sig));
% y(y<0.1)=0;
% y(y>=0.1)=1;
% tmp_sig=tmp_sig.*y;

%tmp_sig=real(bat_signal());
%%%%%%%%%%%%%%%
M=80;
L=2;
[Inew,Iorient]=DTFD(tmp_sig,2,30,M);
figure;tfsapl(tmp_sig,Inew);
Inew=Inew/sum(sum((Inew)));

%tfd_measure(I)