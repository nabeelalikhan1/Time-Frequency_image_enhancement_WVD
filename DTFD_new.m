function [ Inew,Iorient ] = DTFD_new(tmp_sig,alpha1,alpha2,M )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
I= quadtfd(tmp_sig, length(tmp_sig)-1, 1, 'wvd',length(tmp_sig));
[X,Y]=meshgrid(-1:2/M:1,-1:2/M:1);
A=exp((-1/2)*(((25*X).^2)+(25*Y).^2));

I=filter2(A,I);
Iorient=orient_img_calc(abs(I),15,15,15);
%Iorient=orient_img_calc((I),5,5,5);
[x,y]=size(I);
I1=zeros(x+M,y+M);
I1(M/2+1:end-M/2,M/2+1:end-M/2)=I;
% (max(max(Iorient))/pi-.5)*16
% (min(min(Iorient))/pi-.5)*16
Mr=128;
for k=0:Mr
    
    X1=X*cos(pi/2+pi*k/Mr)-Y*sin(pi/2+pi*k/Mr);
    Y1=X*sin(pi/2+pi*k/Mr)+Y*cos(pi/2+pi*k/Mr);
    
    A(:,:,k+1)=exp((-1/2)*(((alpha1*X1).^2)+(alpha2*Y1).^2));
    
end

for i=M/2+1:x+M/2
    for j=M/2+1:y+M/2
        a=round((Iorient(i-M/2,j-M/2)/pi-0.5)*Mr)+1;
%         a
        B=I1(i-M/2:i+M/2,j-M/2:j+M/2);
         
        Inew(i-M/2,j-M/2)=sum(sum(A(:,:,a).*B));
    end
end
end

