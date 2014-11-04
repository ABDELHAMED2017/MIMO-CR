NrxBS = 4;
NtxBS = 4;
NrxD2D = 4;
NtxD2D = 4;

P = 100;
Hba = randn(NrxD2D,NtxD2D);
He = randn(NrxBS,NtxBS);
Hb = randn(NrxBS,NtxD2D);

B = eye(NrxBS);
C = P/NtxBS*eye(NtxBS);

i = 0;
rate_range = 0.1:2:40;
for r_min = rate_range(:)'
    i = i+1;
    Bmax = max_B(B,He,P,r_min,'sum-power');
    [rate_dl_sp,R,Q] = min_max_dl(He,C,Bmax,'sum-power');
    [rate_ul_sp,Om,Si] = min_max_ul(He,C,Bmax,'sum-power');
    rate_al_sp = log(det(eye(NrxBS) + trace(C)/trace(Bmax)*He*He'));
    [rate_dl_sp,rate_ul_sp,rate_al_sp]
    if trace(Bmax) >= NrxBS
        rate_sc_sp(i) = max_scnd(Hba,Hb,P,Bmax,'sum-power');
    else
        rate_sc_sp(i) = 0;
    end
    

    Bmax = max_B(B,He,P,r_min,'per-antenna');
    [rate_dl_pa,R,Q] = min_max_dl(He,C,Bmax,'per-antenna');
    [rate_ul_pa,Om,Si] = min_max_ul(He,C,Bmax,'per-antenna');
    [rate_dl_pa,rate_dl_pa]
    if Bmax(1,1) >= 1
        rate_sc_pa(i) = max_scnd(Hba,Hb,P,Bmax,'per-antenna');
    else
        rate_sc_pa(i) = 0;
    end
    

    Bmax = max_B(B,He,P,r_min,'shape');
    [rate_dl_sh,R,Q] = min_max_dl(He,C,Bmax,'shape');
    [rate_ul_sh,Om,Si] = min_max_ul(He,C,Bmax,'shape');
    [rate_dl_sh,rate_dl_sh]
    if Bmax(1,1) >= 1
        rate_sc_sh(i) = max_scnd(Hba,Hb,P,Bmax,'shape');
    else
        rate_sc_sh(i) = 0;
    end    
end

plot(rate_range,rate_sc_sp,'blue')
hold on
plot(rate_range,rate_sc_pa,'green')
plot(rate_range,rate_sc_sh,'red')
ylabel('Secondary Network')
xlabel('Primary Network')
legend('Sum-Power','Per Antenna','Shaping')
