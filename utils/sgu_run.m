function [X_sim, Y_sim] = sgu_run(modelname,order_app,shocks)
rbc_sv_params;

% The DSGE model
[f,x,xp,y,yp,eta] = feval(modelname,STDA,STDSIG);

% Compute analytical derivatives up to 3 order if needed
[fx,fxp,fy,fyp,...                                                          %first order derivatives
    fxx,fxxp,fxy,fxyp, fxpxp,fxpy,fxpyp, fyy,fyyp, fypyp,...                %second order derivatives
    fxxx,fxxxp,fxxy,fxxyp , fxxpxp,fxxpy,fxxpyp, fxyy,fxyyp, fxypyp, ...    %third order derivatives related to: fxx,fxxp,fxy,fxyp   
    fxpxpxp,fxpxpy,fxpxpyp, fxpyy,fxpyyp, fxpypyp, ...                      %third order derivatives related to: fxpxp,fxpy,fxpyp
    fyyy,fyyyp, fyypyp, ...                                                 %third order derivatives related to: fyy,fyyp
    fypypyp...                                                              %third order derivatives related to: fypyp
    ] = Anal_derivatives(f,x,xp,y,yp,order_app);

% dimensions
n = size(f,1);
nx = size(fx,2);
ny = size(fy,2);
ne = size(shocks,1);

% compute steady-state
Kss = ((1/BETTA-1+DELT)/(Ass*ALPH))^(1/(ALPH-1));
Css = Ass*Kss^ALPH - DELT*Kss;
ln_c_ba = log(Css); ln_c_bap = log(Css); ln_c_cu = log(Css); ln_c_cup = log(Css);
ln_k_ba = log(Kss); ln_k_bap = log(Kss); ln_k_cu = log(Kss); ln_k_cup = log(Kss);
ln_a_ba = log(Ass); ln_a_bap = log(Ass); ln_a_cu = log(Ass); ln_a_cup = log(Ass); 
ln_v_ba = log(1); ln_v_bap = log(1); ln_v_cu = log(1); ln_v_cup = log(1);
ln_siga_ba = log(SIGAss); ln_siga_bap = log(SIGAss); ln_siga_cu = log(SIGAss); ln_siga_cup = log(SIGAss);
eps_a_cu = 0; eps_a_cup = 0;
eps_sig_cu = 0; eps_sig_cup = 0;


%Obtain numerical derivatives of f
num_eval_3rd;
clear f*

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
[SIGy,SIGx]=mom(gx,hx,eta*eta');

%Second-order approximation
if order_app > 1
    [gxx,hxx] = gxx_hxx_noloop(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
else
    gxx = zeros(ny,nx,nx);
    gss = zeros(ny,1);
    hxx = zeros(nx,nx,nx);
    hss = zeros(nx,1);
end

%Third-order approximation
vectorMom3 = zeros(ne,1);
if order_app == 3
    [gxxx,hxxx,gssx,hssx,gsss,hsss] = g_h_3rd(nfx,nfxp,nfy,nfyp,... 
    nfxx,nfxxp,nfxy,nfxyp,nfxpx,nfxpxp,nfxpy,nfxpyp,...
    nfyx,nfyxp,nfyy,nfyyp,nfypx,nfypxp,nfypy,nfypyp,...
    nfxxx,nfxxxp,nfxxy,nfxxyp,nfxxpxp,nfxxpy,nfxxpyp,...
    nfxyy,nfxyyp,nfxypyp,nfxpxpxp,nfxpxpy,nfxpxpyp,...
    nfxpyy,nfxpyyp,nfxpypyp,...
    nfyyy,nfyyyp,nfyypyp,nfypypyp,gx,gxx,gss,hx,hxx,hss,eta,vectorMom3);
else
    gxxx = zeros(ny,nx,nx,nx);
    gssx = zeros(ny,nx);
    gsss = zeros(ny,1);
    hxxx = zeros(nx,nx,nx,nx);    
    hssx = zeros(nx,nx);
    hsss = zeros(nx,1);    
end

%Simulating
sig = 1; % perturbation parameter
[ny,nx] = size(gx);
ne      = size(eta,2);
num_sim = size(shocks,2);

% Allocating memory
Y_sim   = zeros(ny,num_sim);
X_sim   = zeros(nx,num_sim);

% The simulation
xt    = zeros(nx,1);    %The first order terms
AA    = kron(xt,xt);

% Defining matrices
HHxxtil  = 1/2*reshape(hxx,nx,nx^2);
HHxxxtil = 1/6*reshape(hxxx,nx,nx^3);
GGxxtil  = 1/2*reshape(gxx,ny,nx^2);   
GGxxxtil = 1/6*reshape(gxxx,ny,nx^3);   

if order_app == 1
   for t=1:num_sim
      % The state variables
      xt_p  = hx*xt + sig*eta*shocks(:,t);
      X_sim(:,t) = xt_p;

      % The controls
      Y_sim(:,t) = gx*xt_p;
    
      % Updating xt
      xt  = xt_p;
   end
elseif order_app == 2
   for t=1:num_sim
      % The state variables
      xt_p  = hx*xt + sig*eta*shocks(:,t)...
              + HHxxtil*AA + 1/2*sig^2*hss;
      X_sim(:,t) = xt_p;

      % The controls
      AA = kron(xt_p,xt_p);
      Y_sim(:,t) = gx*xt_p + GGxxtil*AA + 1/2*sig^2*gss;
    
      % Updating xt
      xt  = xt_p;
   end
elseif order_app == 3
   for t=1:num_sim
      % The state variables
      xt_p  = hx*xt + sig*eta*shocks(:,t)...
             + HHxxtil*AA + 1/2*sig^2*hss...
             + HHxxxtil*kron(xt,AA)+ 3/6*hssx*sig^2*xt + 1/6*sig^3*hsss;
      X_sim(:,t) = xt_p;    
      
      % The observables
      AA    = kron(xt_p,xt_p);
      Y_sim(:,t) = gx*xt_p + GGxxtil*AA + 1/2*sig^2*gss + ...
                   GGxxxtil*kron(xt_p,AA) + 3/6*gssx*sig^2*xt_p + 1/6*sig^3*gsss;
        
      % Updating xt
      xt= xt_p;
   end
end
