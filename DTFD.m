function [ Inew,Iorient ] = DTFD(tmp_sig,alpha1,alpha2,M )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
I= quadtfd(tmp_sig, length(tmp_sig)-1, 1, 'wvd',length(tmp_sig));
[X,Y]=meshgrid(-1:2/M:1,-1:2/M:1);
A=exp((-1/2)*(((33*X).^2)+(33*Y).^2));
%I=filter2(A,I);

Iorient=orient_img_calc(I,5,5,5);
%Iorient=orient_img_calc_new(I,20);



[x,y]=size(I);
I1=zeros(x+M,y+M);
I1(M/2+1:end-M/2,M/2+1:end-M/2)=I;
for i=M/2+1:length(tmp_sig)+M/2
    for j=M/2+1:length(tmp_sig)+M/2
        X1=X*cos(Iorient(i-M/2,j-M/2))-Y*sin(Iorient(i-M/2,j-M/2));
        Y1=X*sin(Iorient(i-M/2,j-M/2))+Y*cos(Iorient(i-M/2,j-M/2));
        %        X1=X;
        %        Y1=Y;
        A=exp((-1/2)*(((alpha1*X1).^2)+(alpha2*Y1).^2));
        %   A=A.*(1-alpha2*alpha2*Y1.^2);

        
        
        B=I1(i-M/2:i+M/2,j-M/2:j+M/2);
        
        
        Inew(i-M/2,j-M/2)=sum(sum(A.*B));
        
        
        
    end
end
%Inew(Inew<0)=0;
%mesh(Inew)

%Inew=filter2(ones(3,3),Inew);
end

