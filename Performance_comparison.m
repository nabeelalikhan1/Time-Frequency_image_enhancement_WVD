clear all;
close all;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 1:256;
%t = -127:128;
%load seizure;

tmp_sig=4*exp(-0.001*(t-192).^2).*(cos((2*pi*0.07*t)))+2*cos(2*pi*(0.05*t+0.0000015*t.^3))+2*cos(2*pi*(0.075*t+0.0000015*t.^3))+4*exp(-0.001*(t-64).^2).*(cos((2*pi*0.45*t)))+4*exp(-0.001*(t-192).^2).*(cos((2*pi*0.45*t)));%+1*cos(2*pi*(0.45*t));
%tmp_sig =sez_dat(63,:);
%tmp_sig=2*cos(2*pi*(0.05*t+0.0005*t.^2))+2*cos(2*pi*(0.25*t-0.0005*t.^2))+2*exp(-0.001*(t-64).^2).*(cos((2*pi*0.45*t)))+2*exp(-0.001*(t-192).^2).*(cos((2*pi*0.45*t)));%+1*cos(2*pi*(0.45*t));

%tmp_sig=2*cos(2*pi*(0.05*t+0.00075*t.^2))+2*cos(2*pi*(0.15*t+0.00075*t.^2));

%tmp_sig=2*cos(2*pi*(0.05*t+0.0000015*t.^3))+2*cos(2*pi*(0.15*t+0.0000015*t.^3));%+1*cos(2*pi*(0.45*t));

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
%Inew(Inew<0)=0;
(sum(sum(abs(Inew).^(1/L))))^L;
%renyi(Inew,[1:256],[1:256]',3)
%display_orientation(Inew,Iorient,11);

M=60;
I1=HTFD_new1(tmp_sig,3,6,M);    
figure;tfsapl(tmp_sig,I1);
I1=I1/sum(sum((I1)));
(sum(sum(abs(I1).^(1/L))))^L

%renyi(I1,[1:256],[1:256]',3)


I2=HTFD_new1(tmp_sig,2,22,M);
figure;tfsapl(tmp_sig,I2);
I2=I2/sum(sum((I2)));
(sum(sum(abs(I2).^(1/L))))^L

%renyi(I2,[1:256],[1:256]',3)



adaptive_optimal_tfd;
figure;tfsapl(tmp_sig,I_max_new);
I_max_new=I_max_new/sum(sum((I_max_new)));
(sum(sum(abs(I_max_new).^(1/L))))^L

%renyi(I_max_new,[1:256],[1:256]',3)


tfdS=S_method(tmp_sig,10,81);


figure;tfsapl(tmp_sig,tfdS);
tfdS=tfdS/sum(sum((tfdS)));
(sum(sum(abs(tfdS).^(1/L))))^L


%renyi(tfdS,[1:256],[1:256]',3)


tfd_spec=quadtfd(tmp_sig,length(tmp_sig)-1,1,'specx',85,'hamm',length(tmp_sig));
tfd_spec=tfd_spec/sum(sum((tfd_spec)));
(sum(sum(abs(tfd_spec).^(1/L))))^L

%renyi(tfd_spec,[1:256],[1:256]',3)

% tfd_measure(Inew)
%tfd_measure(I)