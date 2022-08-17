% By Martin M. Andreasen, June 8 2009
% This function computes analytical derivatives of the DSGE model
% We use all the symmetry when computing these derivatives. 

function [fx,fxp,fy,fyp,...                                                 %first order derivatives
    fxx,fxxp,fxy,fxyp, fxpxp,fxpy,fxpyp, fyy,fyyp, fypyp,...                %second order derivatives
    fxxx,fxxxp,fxxy,fxxyp , fxxpxp,fxxpy,fxxpyp, fxyy,fxyyp, fxypyp, ...    %third order derivatives related to: fxx,fxxp,fxy,fxyp   
    fxpxpxp,fxpxpy,fxpxpyp, fxpyy,fxpyyp, fxpypyp, ...                      %third order derivatives related to: fxpxp,fxpy,fxpyp
    fyyy,fyyyp, fyypyp, ...                                                 %third order derivatives related to: fyy,fyyp
    fypypyp...                                                              %third order derivatives related to: fypyp
    ] = Anal_derivatives(f,x,xp,y,yp,order_app);

% some indices
nx  = size(x,2);
ny  = size(y,2);
n   = size(f,1);

% The first order derivatives
fx  = jacobian(f,x);
fxp = jacobian(f,xp);
fy  = jacobian(f,y);
fyp = jacobian(f,yp);

% The second order derivatives if needed
if order_app > 1
    % Related to fx
    fxx   = reshape(jacobian(fx(:),x) ,n,nx,nx);
    fxxp  = reshape(jacobian(fx(:),xp),n,nx,nx);    
    fxy   = reshape(jacobian(fx(:),y) ,n,nx,ny);    
    fxyp  = reshape(jacobian(fx(:),yp),n,nx,ny);    
    % Related to fxp
    fxpxp = reshape(jacobian(fxp(:),xp),n,nx,nx);    
    fxpy  = reshape(jacobian(fxp(:),y) ,n,nx,ny);    
    fxpyp = reshape(jacobian(fxp(:),yp),n,nx,ny);    
    % Related to fy
    fyy   = reshape(jacobian(fy(:),y) ,n,ny,ny);    
    fyyp  = reshape(jacobian(fy(:),yp),n,ny,ny);
    % Related to fyp    
    fypyp = reshape(jacobian(fyp(:),yp),n,ny,ny);
else
    %If second order approximation is not needed, we set all second order
    %derivatives to zero
    % Related to fx
    fxx   = zeros(n,nx,nx);
    fxxp  = zeros(n,nx,nx);    
    fxy   = zeros(n,nx,ny);    
    fxyp  = zeros(n,nx,ny);    
    % Related to fxp   
    fxpxp = zeros(n,nx,nx);    
    fxpy  = zeros(n,nx,ny);    
    fxpyp = zeros(n,nx,ny);    
    % Related to fy    
    fyy   = zeros(n,ny,ny);    
    fyyp  = zeros(n,ny,ny);
    % Related to fyp    
    fypyp = zeros(n,ny,ny);    
end

% The third order derivatives if needed
if order_app > 2
    % Related to fxx
    fxxx   = reshape(jacobian(fxx(:),x) ,n,nx,nx,nx);   
    fxxxp  = reshape(jacobian(fxx(:),xp),n,nx,nx,nx);       
    fxxy   = reshape(jacobian(fxx(:),y) ,n,nx,nx,ny);           
    fxxyp  = reshape(jacobian(fxx(:),yp),n,nx,nx,ny);               
    
    % Related to fxxp
    fxxpxp = reshape(jacobian(fxxp(:),xp),n,nx,nx,nx); 
    fxxpy  = reshape(jacobian(fxxp(:),y) ,n,nx,nx,ny);     
    fxxpyp = reshape(jacobian(fxxp(:),yp),n,nx,nx,ny);         
    
    % Related to fxy    
    fxyy   = reshape(jacobian(fxy(:),y) ,n,nx,ny,ny);
    fxyyp  = reshape(jacobian(fxy(:),yp),n,nx,ny,ny);    

    % Related to fxyp
    fxypyp = reshape(jacobian(fxyp(:),yp),n,nx,ny,ny);    
    
    % Related to fxpxp
    fxpxpxp= reshape(jacobian(fxpxp(:),xp),n,nx,nx,nx);   
    fxpxpy = reshape(jacobian(fxpxp(:),y) ,n,nx,nx,ny);       
    fxpxpyp= reshape(jacobian(fxpxp(:),yp),n,nx,nx,ny);
    
    % Related to fxpy
    fxpyy  = reshape(jacobian(fxpy(:),y) ,n,nx,ny,ny);       
    fxpyyp = reshape(jacobian(fxpy(:),yp),n,nx,ny,ny);           
    
    % Related to fxpyp
    fxpypyp= reshape(jacobian(fxpyp(:),yp),n,nx,ny,ny);       
    
    % Related to fyy
    fyyy   = reshape(jacobian(fyy(:),y) ,n,ny,ny,ny);       
    fyyyp  = reshape(jacobian(fyy(:),yp),n,ny,ny,ny);           
    
    % Related to fyyp
    fyypyp = reshape(jacobian(fyyp(:),yp),n,ny,ny,ny);       
    
    % Related to fypyp
    fypypyp= reshape(jacobian(fypyp(:),yp),n,ny,ny,ny);       
else
    %If third order approximation is not needed, we set all third order
    %derivatives to zero
    % Related to fxx
    fxxx   = zeros(n,nx,nx,nx);   
    fxxxp  = zeros(n,nx,nx,nx);       
    fxxy   = zeros(n,nx,nx,ny);           
    fxxyp  = zeros(n,nx,nx,ny);               
    % Related to fxxp
    fxxpxp = zeros(n,nx,nx,nx); 
    fxxpy  = zeros(n,nx,nx,ny);     
    fxxpyp = zeros(n,nx,nx,ny);         
    % Related to fxy    
    fxyy   = zeros(n,nx,ny,ny);
    fxyyp  = zeros(n,nx,ny,ny);    
    % Related to fxyp
    fxypyp = zeros(n,nx,ny,ny);    
    % Related to fxpxp
    fxpxpxp= zeros(n,nx,nx,nx);   
    fxpxpy = zeros(n,nx,nx,ny);       
    fxpxpyp= zeros(n,nx,nx,ny);
    % Related to fxpy
    fxpyy  = zeros(n,nx,ny,ny);       
    fxpyyp = zeros(n,nx,ny,ny);           
    % Related to fxpyp
    fxpypyp= zeros(n,nx,ny,ny);       
    % Related to fyy
    fyyy   = zeros(n,ny,ny,ny);       
    fyyyp  = zeros(n,ny,ny,ny);           
    % Related to fyyp
    fyypyp = zeros(n,ny,ny,ny);       
    % Related to fypyp
    fypypyp= zeros(n,ny,ny,ny);       
end










    