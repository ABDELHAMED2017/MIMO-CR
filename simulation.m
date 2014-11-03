Nrx = 4;
Ntx = 4;
P = 10;
Hba = randn(Nrx,Ntx);
He = randn(Nrx,Ntx);
Hb = randn(Nrx,Ntx);

B = eye(Nrx);
C = P/Ntx*eye(Ntx);
c = norm(C,'fro');
i = 0;
rate_range = 0.1:2:40;
for r_min = rate_range(:)'
    i = i+1;
    Bmax = max_B(B,Hb,P,r_min,'sum-power');
    [rate_dl_sp,R,Q] = min_max_dl(Hb,C,Bmax,'sum-power');
    [rate_ul_sp,Om,Si] = min_max_ul(Hb,C,Bmax,'sum-power');
    rate_al_sp = log(det(eye(Nrx) + trace(C)/trace(Bmax)*Hb*Hb'));
    [rate_dl_sp,rate_ul_sp,rate_al_sp]
    if trace(Bmax) >= Nrx
        rate_sc_sp(i) = max_scnd(Hba,Hb,P,Bmax,'sum-power');
    else
        rate_sc_sp(i) = 0;
    end
    

    Bmax = max_B(B,Hb,P,r_min,'per-antenna');
    [rate_dl_pa,R,Q] = min_max_dl(Hb,C,Bmax,'per-antenna');
    [rate_ul_pa,Om,Si] = min_max_ul(Hb,C,Bmax,'per-antenna');
    [rate_dl_pa,rate_dl_pa]
    rate_sc_pa(i) = max_scnd(Hba,Hb,P,Bmax,'per-antenna');
    if Bmax(1,1) >= 1
        rate_sc_pa(i) = max_scnd(Hba,Hb,P,Bmax,'per-antenna');
    else
        rate_sc_pa(i) = 0;
    end
    

    Bmax = max_B(B,Hb,P,r_min,'shape');
    [rate_dl_sh,R,Q] = min_max_dl(Hb,C,Bmax,'shape');
    [rate_ul_sh,Om,Si] = min_max_ul(Hb,C,Bmax,'shape');
    [rate_dl_sh,rate_dl_sh]
    rate_sc_sh(i) = max_scnd(Hba,Hb,P,Bmax,'shape');
    if Bmax(1,1) >= 1
        rate_sc_sh(i) = max_scnd(Hba,Hb,P,Bmax,'shape');
    else
        rate_sc_sh(i) = 0;
    end    
end

plot(rate_range,rate_sc_pa,'blue')
hold on
plot(rate_range,rate_sc_pa,'green')
plot(rate_range,rate_sc_sh,'red')
ylabel('Secondary Network')
xlabel('Primary Network')
legend('Sum-Power','Per Antenna','Shaping')
