function [f,x,xp,y,yp,eta] = rbc_sv_FL(STDA,STDSIG)
syms ALPH BETTA GAMA DELT RHOA RHOSIG Ass SIGAss
syms ln_c_ba  ln_k_ba   ln_a_ba   ln_v_ba   ln_siga_ba
syms ln_c_bap ln_k_bap  ln_a_bap  ln_v_bap  ln_siga_bap
syms ln_c_cu  ln_k_cu   ln_a_cu   ln_v_cu   ln_siga_cu
syms ln_c_cup ln_k_cup  ln_a_cup  ln_v_cup  ln_siga_cup
syms eps_a_cu  eps_sig_cu
syms eps_a_cup eps_sig_cup

fc = -exp(ln_c_cu)^(-GAMA) + BETTA*(exp(ln_c_cup))^(-GAMA)*(exp(ln_a_cup)*ALPH*exp(ln_k_cup)^(ALPH-1) + 1 - DELT );
fk = -exp(ln_c_cu) -exp(ln_k_cup) + exp(ln_a_cu)*exp(ln_k_cu)^ALPH + (1-DELT)*exp(ln_k_cu);
fa = -(ln_a_cu - log(Ass)) + RHOA*(ln_a_ba - log(Ass)) + exp(ln_siga_cu)*eps_a_cu;
fsiga = -(ln_siga_cup - log(SIGAss)) + RHOSIG*(ln_siga_cu - log(SIGAss));
feps_a = -eps_a_cup;
fa_ba = -ln_a_bap + ln_a_cu;
f = [fc;fk;fa;fsiga;feps_a;fa_ba];

y  = [ln_c_cu  ln_a_cu];
yp = [ln_c_cup ln_a_cup]; 
x  = [ln_k_cu  ln_a_ba  ln_siga_cu  eps_a_cu];
xp = [ln_k_cup ln_a_bap ln_siga_cup eps_a_cup];

eta = zeros(size(x,1),2);
eta(3,1) = STDSIG;
eta(4,2) = STDA;