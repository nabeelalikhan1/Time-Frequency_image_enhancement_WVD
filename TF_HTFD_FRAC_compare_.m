clear all;
close all;
dis=0;
lll=1;
LL=3;
N_S=20;
type=3;
[t,x,fs,IF_O]=signal_type_new(type);
N_C=3;
if type==4
    N_C=2;
end
xn=x;
%type=15
for snr=0:3:12
    
    
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        I_max_new=tfr_stft_high(hilbert(x'));
        % figure; imshow(1-I_max_new,[]);
        %                 figure;imshow(1-I_fixed_wind,[]);
        %IF_COMPUTE_new_edge_link;
        IF_COMPUTE;
        
    end
    var_snr_AFS(lll,:)=(mean(mse1))';
%     
%     
%     
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        
        if type==3
        I_max_new = quadtfd( x, length(x)-1, 1, 'specx', 73, 'hamm',256);
        elseif type==4
        I_max_new = quadtfd( x, length(x)-1, 1, 'specx', 91, 'hamm',256);
        end
            I_max_new=imresize(I_max_new,[length(x) length(x)]);
        
        I_max_new(I_max_new<0)=0;
        
        %IF_COMPUTE_new_edge_link;
        IF_COMPUTE;
        
        
    end
    var_snr_spec(lll,:)=(mean(mse1))';
    
    
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        wvdz=wvd(x,length(x)-1,1,2^nextpow2(length(x)));
        if type==3
        g=extnd_mbd(0.05,0.125,0.5,length(x));
        elseif type==4
             g=extnd_mbd(0.055,0.125,1,length(x));
        end
                gsig=ifft(fft(wvdz.').');

        smg=gsig.*g;
        I_max_new=real(fft(ifft(smg.').'));
        I_max_new=imresize(I_max_new,[length(x) length(x)]);
        I_max_new(I_max_new<0)=0;
        %IF_COMPUTE_new_edge_link;
        IF_COMPUTE;
        
        
    end
    var_snr_embd(lll,:)=(mean(mse1))';
    
    
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        wvdz=wvd(x,length(x)-1,1,2^nextpow2(length(x)));
        if type==3
        g=extnd_mbd(0.03,1,0.5,length(x));
        elseif type==4
             g=extnd_mbd(0.05,1,1,length(x));
        end
                gsig=ifft(fft(wvdz.').');

        smg=gsig.*g;
        I_max_new=real(fft(ifft(smg.').'));
%         [I_max_new Iorient]=DTFD(x,2,22,60);
        I_max_new=imresize(I_max_new,[length(x) length(x)]);
        I_max_new(I_max_new<0)=0;
        %IF_COMPUTE_new_edge_link;
        IF_COMPUTE;
        
        
    end
    var_snr_mbd(lll,:)=(mean(mse1))';
    
    
    
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        
        
        %I_max_new = HTFD_new1(x, 2,30*2,85);%2 20
     %   I_max_new = HTFD_new1(x, 2,20,85);%2 20
                [I_max_new Iorient1]=DTFD(x,2,22,80);

        I_max_new(I_max_new<0)=0;
        I_max_new=imresize(I_max_new,[length(x) length(x)]);
        
        %IF_COMPUTE_new_edge_link;
        IF_COMPUTE;
        
    end
    var_snr_HTFD(lll,:)=(mean(mse1))';
%     
    for k1=1:N_S
        [t,x,fs,IF_O]=signal_type_new(type);
        
        x=awgn(x,snr,'measured',lll*200+k1+3+10);
        [ am I]=wvd1(x',length(x'));
        
        %%%% 150
        C=5;
        D=0.15;
        E=0.075;
        
        %%% 173
        %D=0.12;
%         C=1;
%         D=0.1;
%         E=0.0251;
        g=cskabedbelchourini(length(x),C,D, E);
        am=am.*g;
        I_max_new = (fft(ifftshift(am,1), [], 1));
        I_max_new=  ifft(fftshift(I_max_new,2), [], 2);
        I_max_new=real(I_max_new);
        I_max_new(I_max_new<0)=0;
        
        
        %IF_COMPUTE_new_edge_link;
        IF_COMPUTE;
        
    end
    var_snr_CSK(lll,:)=(mean(mse1))';
    
    lll=lll+1;
end



% snr=-5:3:15;
% 
% display('frac tfd');
% log10(var_snr_AFS)
% save new_mse_AFS_type2.mat var_snr_AFS
% 
% 
% display('MBD');
% log10(var_snr_embd)
% save new_mse_mbd_type2.mat var_snr_mbd
% 
% display('HTFD');
% log10(var_snr_HTFD)
% save new_mse_mbd_type2.mat var_snr_HTFD


snr=0:3:12;
for i=1:1
%      figure;
%     plot(snr,mean(log10(var_snr_mbd')),'-co','linewidth',3);
% %    
    hold on;
    plot(snr,mean(log10(var_snr_spec')),':gs','linewidth',3);
    
%     hold on;
%     plot(snr,mean(log10(var_snr_embd')),'-.b+','linewidth',3);
    hold on;
    plot(snr,mean(log10(var_snr_CSK')),'--kd','linewidth',3);
    
    hold on;
    plot(snr, mean(log10(var_snr_AFS')),'-rh','linewidth',3);

    hold on;
    plot(snr, mean(log10(var_snr_HTFD')),'-.b+','linewidth',3);

    xlabel('Signal to Noise Ratio');
    ylabel('log10(Mean Square Error)');
    axis([min(snr) max(snr)  -5  0])
end


for i=1:3
     figure;
%     plot(snr,(log10(var_snr_mbd(:,i))),'-co','linewidth',3);
%    
   % hold on;
    plot(snr,(log10(var_snr_spec(:,i))),':gs','linewidth',3);
    
%     hold on;
%     plot(snr,(log10(var_snr_embd(:,i))),'-.b+','linewidth',3);
    hold on;
    plot(snr,(log10(var_snr_CSK(:,i))),'--kd','linewidth',3);
    
    hold on;
    plot(snr, (log10(var_snr_AFS(:,i))),'-rh','linewidth',3);

    hold on;
    plot(snr, (log10(var_snr_HTFD(:,i))),'-.b+','linewidth',3);

    xlabel('Signal to Noise Ratio');
    ylabel('log10(Mean Square Error)');
    axis([min(snr) max(snr)  -5  0])
end
