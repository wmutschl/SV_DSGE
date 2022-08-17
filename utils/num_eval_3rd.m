% By Martin M. Andreasen, June 8 2009
% This function evalutes the model and all non-symmetric derivatives up to 
% third order in steady state. 
% This file is a modified version of the code of the same name made by
% Stephanie Schmitt-Grohe and Martin Uribe

% The model evaluated in the steady state
nf     = zeros(n,1);
nf(:)  = eval(f(:));

% The first order derivatives of the model
nfx    = zeros(n,nx);
nfx(:) = eval(fx(:));
nfxp   = zeros(n,nx);
nfxp(:)= eval(fxp(:));
nfy    = zeros(n,ny);
nfy(:) = eval(fy(:));
nfyp   = zeros(n,ny);
nfyp(:)= eval(fyp(:));

% The second order derivatives of the model
if order_app > 1
    % Derivatives related to fx
    nfxx     = zeros(n,nx,nx);
    nfxx(:)  = eval(fxx(:));
    nfxxp    = zeros(n,nx,nx);
    nfxxp(:) = eval(fxxp(:));
    nfxy     = zeros(n,nx,ny);
    nfxy(:)  = eval(fxy(:));
    nfxyp    = zeros(n,nx,ny);
    nfxyp(:) = eval(fxyp(:));    
   
    % Derivatives related to fxp
    % Using Young's theorem
    nfxpx    = zeros(n,nx,nx);
    for i=1:n
        nfxpx(i,:,:) = squeeze(nfxxp(i,:,:))';
    end
    nfxpxp   = zeros(n,nx,nx);
    nfxpxp(:)= eval(fxpxp(:));
    nfxpy    = zeros(n,nx,ny);    
    nfxpy(:) = eval(fxpy(:));
    nfxpyp   = zeros(n,nx,ny);    
    nfxpyp(:)= eval(fxpyp(:));    
    
    % Derivatives related to fy
    % Using Young's theorem
    nfyx = zeros(n,ny,nx);
    for i=1:n
        nfyx(i,:,:) = squeeze(nfxy(i,:,:))';
    end
    nfyxp = zeros(n,ny,nx);
    for i=1:n
        nfyxp(i,:,:) = squeeze(nfxpy(i,:,:))';
    end
    nfyy     = zeros(n,ny,ny);
    nfyy(:)  = eval(fyy(:));
    nfyyp    = zeros(n,ny,ny);
    nfyyp(:) = eval(fyyp(:));    
    
    % Derivatives related to fyp
    % Using Young's theorem
    nfypx = zeros(n,ny,nx);
    for i=1:n
        nfypx(i,:,:) = squeeze(nfxyp(i,:,:))';
    end
    nfypxp = zeros(n,ny,nx);
    for i=1:n
        nfypxp(i,:,:) = squeeze(nfxpyp(i,:,:))';
    end
    nfypy = zeros(n,ny,ny);
    for i=1:n
        nfypy(i,:,:) = squeeze(nfyyp(i,:,:))';
    end
    nfypyp    = zeros(n,ny,ny);
    nfypyp(:) = eval(fypyp(:));        
else
    %If second order approximation is not needed, we set all second order
    % Derivatives related to fx
    nfxx     = zeros(n,nx,nx);
    nfxxp    = zeros(n,nx,nx);
    nfxy     = zeros(n,nx,ny);
    nfxyp    = zeros(n,nx,ny);
   
    % Derivatives related to fxp
    nfxpx    = zeros(n,nx,nx);
    nfxpxp   = zeros(n,nx,nx);
    nfxpy    = zeros(n,nx,ny);    
    nfxpyp   = zeros(n,nx,ny);    
    
    % Derivatives related to fy
    nfyx      = zeros(n,ny,nx);
    nfyxp     = zeros(n,ny,nx);
    nfyy      = zeros(n,ny,ny);
    nfyyp     = zeros(n,ny,ny);

    % Derivatives related to fyp
    nfypx     = zeros(n,ny,nx);
    nfypxp    = zeros(n,ny,nx);
    nfypy     = zeros(n,ny,ny);
    nfypyp    = zeros(n,ny,ny);
end 

% The third order derivatives of the model
if order_app > 2
    % Derivatives related to fxx
    nfxxx      = zeros(n,nx,nx,nx);
    nfxxx(:)   = eval(fxxx(:));
    nfxxxp     = zeros(n,nx,nx,nx);
    nfxxxp(:)  = eval(fxxxp(:));    
    nfxxy      = zeros(n,nx,nx,ny);
    nfxxy(:)   = eval(fxxy(:));        
    nfxxyp     = zeros(n,nx,nx,ny);
    nfxxyp(:)  = eval(fxxyp(:));

    % Derivatives related to fxxp
    nfxxpxp    = zeros(n,nx,nx,nx);
    nfxxpxp(:) = eval(fxxpxp(:));
    nfxxpy     = zeros(n,nx,nx,ny);
    nfxxpy(:)  = eval(fxxpy(:));
    nfxxpyp    = zeros(n,nx,nx,ny);
    nfxxpyp(:) = eval(fxxpyp(:));
    
    % Derivatives related to fxy
    nfxyy      = zeros(n,nx,ny,ny);
    nfxyy(:)   = eval(fxyy(:));
    nfxyyp     = zeros(n,nx,ny,ny);
    nfxyyp(:)  = eval(fxyyp(:));
    
    % Derivatives related to fxyp
    nfxypyp    = zeros(n,nx,ny,ny);    
    nfxypyp(:) = eval(fxypyp(:));

    % Related to fxpxp
    nfxpxpxp   = zeros(n,nx,nx,nx);   
    nfxpxpxp(:)= eval(fxpxpxp(:));
    nfxpxpy    = zeros(n,nx,nx,ny);       
    nfxpxpy(:) = eval(fxpxpy(:));
    nfxpxpyp   = zeros(n,nx,nx,ny);
    nfxpxpyp(:)= eval(fxpxpyp(:));
    
    % Related to fxpy
    nfxpyy     = zeros(n,nx,ny,ny);       
    nfxpyy(:)  = eval(fxpyy(:));
    nfxpyyp    = zeros(n,nx,ny,ny);           
    nfxpyyp(:) = eval(fxpyyp(:));
    
    % Related to fxpyp
    nfxpypyp   = zeros(n,nx,ny,ny);       
    nfxpypyp(:)= eval(fxpypyp(:));
    
    % Related to fyy
    nfyyy      = zeros(n,ny,ny,ny);       
    nfyyy(:)   = eval(fyyy(:));
    nfyyyp     = zeros(n,ny,ny,ny);           
    nfyyyp(:)  = eval(fyyyp(:));
    
    % Related to fyyp
    nfyypyp    = zeros(n,ny,ny,ny);       
    nfyypyp(:) = eval(fyypyp(:));           
    
    % Related to fypyp
    nfypypyp   = zeros(n,ny,ny,ny);      
    nfypypyp(:)= eval(fypypyp(:));
else
    %If third order approximation is not needed, we set all third order
    %derivatives to zero
    % Derivatives related to fxx
    nfxxx      = zeros(n,nx,nx,nx);
    nfxxxp     = zeros(n,nx,nx,nx);
    nfxxy      = zeros(n,nx,nx,ny);
    nfxxyp     = zeros(n,nx,nx,ny);

    % Derivatives related to fxxp
    nfxxpxp    = zeros(n,nx,nx,nx);
    nfxxpy     = zeros(n,nx,nx,ny);
    nfxxpyp    = zeros(n,nx,nx,ny);
    
    % Derivatives related to fxy
    nfxyy      = zeros(n,nx,ny,ny);
    nfxyyp     = zeros(n,nx,ny,ny);
    
    % Derivatives related to fxyp
    nfxypyp    = zeros(n,nx,ny,ny);    

    % Related to fxpxp
    nfxpxpxp   = zeros(n,nx,nx,nx);   
    nfxpxpy    = zeros(n,nx,nx,ny);       
    nfxpxpyp   = zeros(n,nx,nx,ny);
    
    % Related to fxpy
    nfxpyy     = zeros(n,nx,ny,ny);       
    nfxpyyp    = zeros(n,nx,ny,ny);           
    
    % Related to fxpyp
    nfxpypyp   = zeros(n,nx,ny,ny);       
    
    % Related to fyy
    nfyyy      = zeros(n,ny,ny,ny);       
    nfyyyp     = zeros(n,ny,ny,ny);           
    
    % Related to fyyp
    nfyypyp    = zeros(n,ny,ny,ny);       
    
    % Related to fypyp
    nfypypyp   = zeros(n,ny,ny,ny);      
end