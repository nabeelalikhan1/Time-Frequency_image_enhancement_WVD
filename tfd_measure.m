function mm_val = tfd_measure(tfd)

% For signals with two LFM components.
% Several points are disregarded at the beginning and end of TFD.

% Written by Shiying Dong, 2012

tfd = tfd./max(tfd(:)); % tfd normalization
tfd = abs(tfd);
[F_data, T_data] = size(tfd);

N_dis = 16; %%% disregard several slices at the beginning and end of tfd
for nn = N_dis+1:T_data-N_dis
    timeslice = tfd(:,nn);
    tmp1 = diff(timeslice);
    tmp_mark = tmp1(2:end).*tmp1(1:end-1);
    tmp_peak = find(tmp_mark<=0 & tmp1(1:end-1)>0)+1;
    tmp_tip = find(tmp_mark<0)+1;
    
%     plot(tfd(:,nn),'.-')
    
    % component IF and IAmp
    tmp_thr = max(timeslice)*0.3;
    tmp1_IF = tmp_peak(timeslice(tmp_peak)>tmp_thr);
    if length(tmp1_IF) == 2
        IF_est(nn,:) = tmp1_IF;
        IAmp_est(nn,:) = timeslice(tmp1_IF);
    else
        IF_est(nn,:) = tmp1_IF([1 length(tmp1_IF)]);
        IAmp_est(nn,:) = timeslice(tmp1_IF([1 length(tmp1_IF)]));
    end
    
    for N_comp = 1:2 %%% number of components
        tmp3_IF = IF_est(nn,N_comp);
        % sidelobe amplitudes
        if isempty(tmp_peak(tmp_peak<tmp3_IF))
            sidelope_L(nn,N_comp) = 0;
            Amp_sidelope_L(nn,N_comp) = 0;
        else
            sidelope_L(nn,N_comp) = max(tmp_peak(tmp_peak<tmp3_IF));
            Amp_sidelope_L(nn,N_comp) = timeslice(sidelope_L(nn,N_comp));
        end
        if isempty(tmp_peak(tmp_peak>tmp3_IF))
            sidelope_R(nn,N_comp) = 0;
            Amp_sidelope_R(nn,N_comp) = 0;
        else
            sidelope_R(nn,N_comp) = min(tmp_peak(tmp_peak>tmp3_IF));
            Amp_sidelope_R(nn,N_comp) = timeslice(sidelope_R(nn,N_comp));
        end
        
        % bandwidth
        if isempty(tmp_tip(tmp_tip<tmp3_IF))
            tmp_st = 1;
        else
            tmp_st = max(tmp_tip(tmp_tip<tmp3_IF));
        end
        if isempty(tmp_tip(tmp_tip>tmp3_IF))
            tmp_en = 128;
        else
            tmp_en = min(tmp_tip(tmp_tip>tmp3_IF));
        end
        comp_envelope = timeslice(tmp_st:tmp_en)-IAmp_est(nn,N_comp)/sqrt(2);
        comp_envelope_interp = interp1(1:length(comp_envelope),comp_envelope,1:0.001:length(comp_envelope),'cubic');
        [interp_peak_val  interp_peak] = max(comp_envelope_interp);
        [bdwdth_L_val  bdwdth_L] = min(abs(comp_envelope_interp(1:interp_peak)));
        [bdwdth_R_val  bdwdth_R] = min(abs(comp_envelope_interp(interp_peak:end)));
        bdwdth(nn,N_comp) = abs(bdwdth_R+interp_peak-1-bdwdth_L)*0.001/F_data/2;  %****
    end
    % cross-term
    tmp_xterm = round(sum(IF_est(nn,:))/2);
    [xterm(nn) xterm_indx(nn)] = max(timeslice(tmp_xterm-2:tmp_xterm+2));
    
    %     test
%     plot(tfd(:,nn),'.-')
%     hold on
%     plot(sidelope_L(nn,:),Amp_sidelope_L(nn,:),'*b')
%     plot(sidelope_R(nn,:),Amp_sidelope_R(nn,:),'*k')
%     plot(IF_est(nn,:),IAmp_est(nn,:),'ro')
%     plot(xterm_indx(nn)+tmp_xterm-2-1,xterm(nn),'^')
%     hold off
%     nn
%     input ''
end
% measure

Amp_sidelope_L(1:N_dis,:) = [];
Amp_sidelope_R(1:N_dis,:) = [];
IAmp_est(1:N_dis,:) = [];
xterm(1:N_dis) = [];
IF_est(1:N_dis,:) = [];
bdwdth(1:N_dis,:) = [];

term1 = (Amp_sidelope_L(:,1)/2+Amp_sidelope_R(:,2)/2)./mean(IAmp_est,2);
term2 = xterm'./mean(IAmp_est,2);
IF_est_val = IF_est/(2*T_data);
term3 = sum(bdwdth,2)./2./abs(diff(IF_est_val')');

mm_val = mean(1-(term1+term2/2+term3)/3);

