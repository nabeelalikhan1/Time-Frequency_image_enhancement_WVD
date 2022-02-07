clear all;
close all;
n=1:256;
sig=zeros(1,256);
sig=2*cos(2*pi*n*0.05);
% sig(30)=sig(30)+20;
% sig(70)=sig(70)+20;
% sig(110)=sig(110)+20;
% sig(150)=sig(150)+20;
% sig(190)=sig(190)+20;
% sig(230)=sig(230)+20;
sigm=2;
B = fir1(16*2,0.08,'high');

s1=10*exp(-((n-30).^2)/sigm);%.*cos(0.1*pi*n);
s1=filter(B,1,s1);

s2=10*exp(-((n-70).^2)/sigm);%.*cos(0.1*pi*n);
s2=filter(B,1,s2);

s3=10*exp(-((n-110).^2)/sigm);%.*cos(0.1*pi*n);
s3=filter(B,1,s3);

s4=10*exp(-((n-150).^2)/sigm);%.*cos(0.1*pi*n);
s4=filter(B,1,s4);

s5=10*exp(-((n-190).^2)/sigm);%.*cos(0.1*pi*n);
s5=filter(B,1,s5);

s6=10*exp(-((n-230).^2)/sigm);%.*cos(0.1*pi*n);
s6=filter(B,1,s6);

s=s1+s2+s3+s4+s5+s6;

I=wvd(s1,length(s1)-1)+wvd(s2,length(s2)-1)+wvd(s3,length(s3)-1)+wvd(s4,length(s4)-1)+wvd(s5,length(s5)-1)+wvd(s6,length(s6)-1);
II=wvd(sig,length(sig)-1)+I;
II=filter2(ones(5,5),II);
%sig=sig+20*s.*cos(0.15*pi*n);
%sig=sig-mean(sig);
figure;tfsapl(sig,II)
sig=sig+s;%.*cos(0.1*pi*n);
%sig=5*s.*cos(0.15*pi*n);
%sig=filter(B,1,sig);
sig=sig-mean(sig);
%I1=HTFD_new1(sig,3,6,81);
I1=DTFD(sig,3,8,65);
figure;
tfsapl(sig,I1)

I=wvd(sig,length(sig)-1);
N=100;
t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^(0.1));
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^(0.2));
h2=h2/sum(h2);
B2=h1'*h2;

I=filter2(B2,I);
figure;
tfsapl(sig,I)

tfd = quadtfd(sig,255,1,'specx',61,'hamm',256);
figure;
tfsapl(sig,tfd)


