% By Martin M. Andreasen, April 26 2010
% This function computes gxxx, hxxx, gssx, hssx, gsss, and hsss
% The output:
%   - gxxx(beta1,alfa1,alfa2,alfa3) where beta1 = 1,2,...,ny and
%     alfa1,alfa2,alfa3 = 1,2,...nx
%   - hxxx(gama,alfa1,alfa2,alfa3) where gama1 = 1,2,...,nx and
%     alfa1,alfa2,alfa3 = 1,2,...nx
%   - gssx a matrix of order ny * nx
%   - hssx a matrix of order nx * nx
%   - gsss a matrix of order ny * 1
%   - hsss a matrix of order nx * 1

function [gxxx,hxxx,gssx,hssx,gsss,hsss] = g_h_3rd(fx,fxp,fy,fyp,...   % First-order terms
    fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,...                                 % Second-order terms
    fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,...
    fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,...                                % Third-order terms
    fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,...
    fxpyy,fxpyyp,fxpypyp,...
    fyyy,fyyyp,fyypyp,fypypyp,gx,gxx,gss,hx,hxx,hss,eta,vectorMom3)

% some integer constants
[ny,nx]  = size(gx);
ne       = size(eta,2);
n        = nx+ny;
nnxI     = n*nx;
nxnxI	 = nx*nx;
nxnyI	 = nx*ny;
nxnxnyI  = nx*nx*ny;
nxnxnxI  = nx*nx*nx;
nnxnxnxI = n*nx*nx*nx;
nxnxnxnyI= nx*nx*nx*ny;
colQ     = (ny+nx)*(nx^2+nx*(nx-1)*(nx-2)/6);

% We start by computing the remaining matrices for the third order derivatives.
% 1) Derivatives related to fxx
%fxxx  : is computed
%fxxxp : is computed
%fxxy  : is computed
%fxxyp : is computed

% 2) Derivatives related to fxxp
fxxpx     = zeros(n,nx,nx,nx);
for i=1:n
    for j=1:nx
        fxxpx(i,j,:,:) = squeeze(fxxxp(i,j,:,:))';
    end
end
%fxxpxp  : is computed
%fxxpy   : is computed
%fxxpyp  : is computed
    
% 3) Derivatives related to fxy
fxyx      = zeros(n,nx,ny,nx);
for i=1:n
    for j=1:nx
        fxyx(i,j,:,:) = squeeze(fxxy(i,j,:,:))';
    end
end
fxyxp     = zeros(n,nx,ny,nx);
for i=1:n
    for j=1:nx
        fxyxp(i,j,:,:) = squeeze(fxxpy(i,j,:,:))';
    end
end    
%fxyy    : is computed
%fxyyp   : is computed

% 4) Derivatives related to fxyp
fxypx     = zeros(n,nx,ny,nx);
for i=1:n
    for j=1:nx
        fxypx(i,j,:,:) = squeeze(fxxyp(i,j,:,:))';
    end
end
fxypxp    = zeros(n,nx,ny,nx);
for i=1:n
    for j=1:nx
        fxypxp(i,j,:,:) = squeeze(fxxpyp(i,j,:,:))';
    end
end    
fxypy    = zeros(n,nx,ny,ny);
for i=1:n
	for j=1:nx
    	fxypy(i,j,:,:) = squeeze(fxyyp(i,j,:,:))';
    end
end    
%fxypyp   : is computed

% 5) Derivatives related to fxpx
fxpxx    = zeros(n,nx,nx,nx);
for i=1:n
    for j=1:nx
        fxpxx(i,:,:,j) = squeeze(fxxpx(i,:,:,j))';
    end
end
fxpxxp   = zeros(n,nx,nx,nx);
for i=1:n
    for j=1:nx
        fxpxxp(i,:,:,j) = squeeze(fxxpxp(i,:,:,j))';
    end
end
fxpxy    = zeros(n,nx,nx,ny);
for i=1:n
    for j=1:ny
        fxpxy(i,:,:,j) = squeeze(fxxpy(i,:,:,j))';
    end
end    
fxpxyp   = zeros(n,nx,nx,ny);
for i=1:n
    for j=1:ny
        fxpxyp(i,:,:,j) = squeeze(fxxpyp(i,:,:,j))';
    end
end        
    
% 6) Derivatives related to fxpxp
fxpxpx    = zeros(n,nx,nx,nx);
for i=1:n
	for j=1:nx
    	fxpxpx(i,:,j,:) = squeeze(fxxpxp(i,:,j,:))';
    end
end    
%fxpxpxp   : is computed
%fxpxpy    : is computed
%fxpxpyp   : is computed

% 7) Derivatives related to fxpy
fxpyx       = zeros(n,nx,ny,nx);
for i=1:n
    for j=1:nx
        fxpyx(i,:,:,j) = squeeze(fxxpy(i,j,:,:));  
    end
end    
fxpyxp      = zeros(n,nx,ny,nx);
for i=1:n
	for j=1:nx
    	fxpyxp(i,j,:,:) = squeeze(fxpxpy(i,j,:,:))';
    end
end
%fxpyy   : is computed
%fxpyyp  : is computed

% 8) Derivatives related to fxpyp
fxpypx      = zeros(n,nx,ny,nx); 
for i=1:n
    for j=1:ny
    	fxpypx(i,:,j,:) = squeeze(fxypxp(i,:,j,:))';
	end
end 
fxpypxp     = zeros(n,nx,ny,nx); 
for i=1:n
	for j=1:nx
    	fxpypxp(i,j,:,:) = squeeze(fxpxpyp(i,j,:,:))';
	end
end
fxpypy      = zeros(n,nx,ny,ny); 
for i=1:n
	for j=1:nx
    	fxpypy(i,j,:,:) = squeeze(fxpyyp(i,j,:,:))';
	end
end     
%fxpyyp : is computed

% 9) Derivatives related to fyx
fyxx        = zeros(n,ny,nx,nx);
for i=1:n
	for j=1:nx
    	fyxx(i,:,j,:) = squeeze(fxxy(i,:,j,:))';
	end
end 
fyxxp       = zeros(n,ny,nx,nx);
for i=1:n
	for j=1:ny
    	fyxxp(i,j,:,:) = squeeze(fxpxy(i,:,:,j))'; %% rettet
	end
end     
fyxy        = zeros(n,ny,nx,ny);
for i=1:n
	for j=1:ny
    	fyxy(i,:,:,j) = squeeze(fxyy(i,:,:,j))';
	end
end
fyxyp       = zeros(n,ny,nx,ny);
for i=1:n
	for j=1:ny
    	fyxyp(i,:,:,j) = squeeze(fxyyp(i,:,:,j))';
	end
end

% 10) Derivatives related to fyxp
fyxpx       = zeros(n,ny,nx,nx);
for i=1:n
	for j=1:nx
    	fyxpx(i,:,j,:) = squeeze(fxxpy(i,:,j,:))';
	end
end
fyxpxp      = zeros(n,ny,nx,nx);
for i=1:n
	for j=1:nx
    	fyxpxp(i,:,j,:) = squeeze(fxpxpy(i,:,j,:))';
	end
end    
fyxpy       = zeros(n,ny,nx,ny);
for i=1:n
	for j=1:ny
    	fyxpy(i,:,:,j) = squeeze(fxpyy(i,:,:,j))';
	end
end        
fyxpyp      = zeros(n,ny,nx,ny);
for i=1:n
	for j=1:ny
    	fyxpyp(i,:,:,j) = squeeze(fxpyyp(i,:,:,j))';
	end
end        

% 11) Derivatives related to fyy
fyyx        = zeros(n,ny,ny,nx);    
for i=1:n
	for j=1:ny
    	fyyx(i,:,j,:) = squeeze(fxyy(i,:,j,:))';
	end
end        
fyyxp       = zeros(n,ny,ny,nx);    
for i=1:n
	for j=1:ny
    	fyyxp(i,:,j,:) = squeeze(fxpyy(i,:,j,:))';
	end
end        
%fyyy  : is computed
%fyyyp : is computed

% 12) Derivatives related to fyyp    
fyypx       = zeros(n,ny,ny,nx);
for i=1:n
	for j=1:nx
    	fyypx(i,:,:,j) = squeeze(fxyyp(i,j,:,:));
	end
end        
fyypxp      = zeros(n,ny,ny,nx);
for i=1:n
	for j=1:nx
    	fyypxp(i,:,:,j) = squeeze(fxpyyp(i,j,:,:));
	end
end        
fyypy      = zeros(n,ny,ny,ny);
for i=1:n
	for j=1:ny
    	fyypy(i,j,:,:) = squeeze(fyyyp(i,j,:,:))';
	end
end   
%fyypyp  : is computed

% 13) Derivatives related to fypx    
fypxx      = zeros(n,ny,nx,nx);
for i=1:n
	for j=1:nx
    	fypxx(i,:,j,:) = squeeze(fxxyp(i,:,j,:))';
	end
end   
fypxxp     = zeros(n,ny,nx,nx);
for i=1:n
	for j=1:ny
    	fypxxp(i,j,:,:) = squeeze(fxxpyp(i,:,:,j));
	end
end   
fypxy     = zeros(n,ny,nx,ny);    
for i=1:n
	for j=1:ny
    	fypxy(i,j,:,:) = squeeze(fxyyp(i,:,:,j));
	end
end   
fypxyp    = zeros(n,ny,nx,ny);
for i=1:n
	for j=1:ny
    	fypxyp(i,:,:,j) = squeeze(fxypyp(i,:,:,j))';
	end
end   

% 14) Derivatives related to fypxp
fypxpx    = zeros(n,ny,nx,nx);    
for i=1:n
	for j=1:ny
    	fypxpx(i,j,:,:) = squeeze(fxxpyp(i,:,:,j))';
	end
end   
fypxpxp   = zeros(n,ny,nx,nx);    
for i=1:n
	for j=1:nx
    	fypxpxp(i,:,j,:) = squeeze(fxpxpyp(i,:,j,:))';
	end
end   
fypxpy    = zeros(n,ny,nx,ny);    
for i=1:n
	for j=1:ny   
    	fypxpy(i,j,:,:) = squeeze(fxpyyp(i,:,:,j));
	end
end   
fypxpyp   = zeros(n,ny,nx,ny);    
for i=1:n
	for j=1:ny
    	fypxpyp(i,:,:,j) = squeeze(fxpypyp(i,:,:,j))';
	end
end   

% 15) Derivatives related to fypy
fypyx     = zeros(n,ny,ny,nx);
for i=1:n
	for j=1:ny
    	fypyx(i,j,:,:) = squeeze(fxyyp(i,:,:,j))';
	end
end   
fypyxp    = zeros(n,ny,ny,nx);
for i=1:n
	for j=1:ny
    	fypyxp(i,j,:,:) = squeeze(fxpyyp(i,:,:,j))';
	end
end   
fypyy     = zeros(n,ny,ny,ny);    
for i=1:n
	for j=1:ny
    	fypyy(i,j,:,:) = squeeze(fyyyp(i,:,:,j));
	end
end   
fypyyp    = zeros(n,ny,ny,ny);    
for i=1:n
	for j=1:ny
    	fypyyp(i,:,:,j) = squeeze(fyypyp(i,:,:,j))';
	end
end   

% 16) Derivatives related to fypyp    
fypypx    = zeros(n,ny,ny,nx);
for i=1:n
	for j=1:nx
    	fypypx(i,:,:,j) = squeeze(fxypyp(i,j,:,:));
	end
end 
fypypxp   = zeros(n,ny,ny,nx);
for i=1:n
	for j=1:nx
    	fypypxp(i,:,:,j) = squeeze(fxpypyp(i,j,:,:));
	end
end 
fypypy    = zeros(n,ny,ny,ny);
for i=1:n
	for j=1:ny
    	fypypy(i,:,:,j) = squeeze(fyypyp(i,j,:,:));
	end
end 
%fypypyp  : is computed



% Allocating memory
a     = zeros(colQ,1);
ass   = zeros(nnxI,1);
Q     = zeros(colQ,nnxnxnxI);
Qr    = zeros(colQ,colQ);
Qss   = zeros(nnxI,nnxI);
gxxx  = zeros(ny,nx,nx,nx);
gssx  = zeros(ny,nx);
hxxx  = zeros(nx,nx,nx,nx);
hssx  = zeros(nx,nx);

% Auxiliary variable
gx_hx = gx*hx;

% The system we solve to get gxxx and hxxx: Q(colQ,nnxnxnxI))*x(nnxnxnxI,1) + a(colQ,1) = 0
% The index m counts the number of unknowns, i.e. the elements in x
m = 0;
for i=1:n
    % An auxiliary variable
    fyp_gx(1,1:nx) = fyp(i,:)*gx(:,1:nx);
    
    for alfa1=1:nx
        for alfa2=1:alfa1
            for alfa3=1:alfa2
                m = m+1;                    
                % Q1 - adding to the constant a(:,1)
                tmp = 0;
                for beta1=1:ny
                    for beta2=1:ny
                        tmp = tmp +(reshape(fypypyp(i,beta1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                                   +reshape(fypypy (i,beta1,beta2,:),1,ny)*gx(:,alfa3)...
                                   +reshape(fypypxp(i,beta1,beta2,:),1,nx)*hx(:,alfa3)...
                                           +fypypx (i,beta1,beta2,alfa3))*gx_hx(beta2,alfa2)*gx_hx(beta1,alfa1)...
                            + fypyp(i,beta1,beta2)*(hx(:,alfa2)'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)*gx_hx(beta1,alfa1)...
                                   +gx(beta2,:)*squeeze(hxx(:,alfa2,alfa3))*gx_hx(beta1,alfa1)...
                                   +gx_hx(beta2,alfa2)*hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa3)...
                                   +gx_hx(beta2,alfa2)*gx(beta1,:)*squeeze(hxx(:,alfa1,alfa3)));
                    end
                end
                a(m,1) = tmp;

                
                % Q2 - adding to the constant a(:,1)                
                %tmp = 0;
                %for beta1=1:ny
                %    for beta2=1:ny
                %        tmp = tmp +(reshape(fypyyp(i,beta1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                %                   +reshape(fypyy (i,beta1,beta2,:),1,ny)*gx(:,alfa3)...
                %                   +reshape(fypyxp(i,beta1,beta2,:),1,nx)*hx(:,alfa3)...
                %                           +fypyx (i,beta1,beta2,alfa3))*gx(beta2,alfa2)*gx_hx(beta1,alfa1)...
                %             +fypy(i,beta1,beta2)*(gxx(beta2,alfa2,alfa3)*gx_hx(beta1,alfa1)...
                %                   +gx(beta2,alfa2)*hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa3)...
                %                   +gx(beta2,alfa2)*gx(beta1,:)*squeeze(hxx(:,alfa1,alfa3)));                        
                %    end
                %end
                tmp = 0;
                for beta1=1:ny
                    tmp = tmp +gx(:,alfa2)'*(reshape(fypyyp(i,beta1,:,:),ny,ny)*gx_hx(:,alfa3)...
                               +reshape(fypyy (i,beta1,:,:),ny,ny)*gx(:,alfa3)...
                               +reshape(fypyxp(i,beta1,:,:),ny,nx)*hx(:,alfa3)...
                               +reshape(fypyx (i,beta1,:,alfa3),ny,1))*gx_hx(beta1,alfa1)...
                         +reshape(fypy(i,beta1,:),1,ny)*(squeeze(gxx(:,alfa2,alfa3))*gx_hx(beta1,alfa1)...
                         +gx(:,alfa2)*hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa3)...
                         +gx(:,alfa2)*gx(beta1,:)*squeeze(hxx(:,alfa1,alfa3)));
                end
                a(m,1) = a(m,1) + tmp;
                
                
                % Q3 - adding to the constant a(:,1)
                %tmp = 0;
                %for beta1=1:ny
                %    for gama2=1:nx
                %        tmp = tmp +(reshape(fypxpyp(i,beta1,gama2,:),1,ny)*gx_hx(:,alfa3)...
                %                   +reshape(fypxpy (i,beta1,gama2,:),1,ny)*gx(:,alfa3)...
                %                   +reshape(fypxpxp(i,beta1,gama2,:),1,nx)*hx(:,alfa3)...
                %                           +fypxpx (i,beta1,gama2,alfa3))*hx(gama2,alfa2)*gx_hx(beta1,alfa1)...
                %             +fypxp(i,beta1,gama2)*(hxx(gama2,alfa2,alfa3)*gx_hx(beta1,alfa1)...
                %                   +hx(gama2,alfa2)*hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa3)...
                %                   +hx(gama2,alfa2)*gx(beta1,:)*squeeze(hxx(:,alfa1,alfa3)));
                %    end
                %end
                tmp = 0;
                for beta1=1:ny
                    tmp = tmp +hx(:,alfa2)'*(reshape(fypxpyp(i,beta1,:,:),nx,ny)*gx_hx(:,alfa3)...
                               +reshape(fypxpy (i,beta1,:,:),nx,ny)*gx(:,alfa3)...
                               +reshape(fypxpxp(i,beta1,:,:),nx,nx)*hx(:,alfa3)...
                               +reshape(fypxpx (i,beta1,:,alfa3),nx,1))*gx_hx(beta1,alfa1)...
                         +reshape(fypxp(i,beta1,:),1,nx)*(squeeze(hxx(:,alfa2,alfa3))*gx_hx(beta1,alfa1)...
                         +hx(:,alfa2)*hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa3)...
                         +hx(:,alfa2)*gx(beta1,:)*squeeze(hxx(:,alfa1,alfa3)));
                end
                a(m,1) = a(m,1) + tmp;
                
                
                % Q4 - adding to the constant a(:,1)                
                tmp = 0;
                for beta1=1:ny
                   tmp = tmp + (reshape(fypxyp(i,beta1,alfa2,:),1,ny)*gx_hx(:,alfa3)...
                               +reshape(fypxy (i,beta1,alfa2,:),1,ny)*gx(:,alfa3)...
                               +reshape(fypxxp(i,beta1,alfa2,:),1,nx)*hx(:,alfa3)...
                                       +fypxx (i,beta1,alfa2,alfa3))*gx_hx(beta1,alfa1)...
                         +fypx(i,beta1,alfa2)*(hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa3)...
                               +gx(beta1,:)*squeeze(hxx(:,alfa1,alfa3)));
                end
                a(m,1) = a(m,1) + tmp;
                
                
                % Q5 - adding to the constant a(:,1) and given loadings for gxxx
                tmp = 0;
                for beta1=1:ny
                    tmp = tmp + (reshape(fypyp(i,beta1,:),1,ny)*gx_hx(:,alfa3)...
                                +reshape(fypy (i,beta1,:),1,ny)*gx(:,alfa3)...
                                +reshape(fypxp(i,beta1,:),1,nx)*hx(:,alfa3)...
                                        +fypx (i,beta1,alfa3))*hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*hx(:,alfa2)...
                          +fyp(i,beta1)*(hx(:,alfa1)'*squeeze(gxx(beta1,:,:))*squeeze(hxx(:,alfa2,alfa3))...
                                +reshape(hxx(:,alfa1,alfa3),1,nx)*squeeze(gxx(beta1,:,:))*hx(:,alfa2));
                end
                a(m,1) = a(m,1) + tmp;
                % Note, we load into the coefficient matrix by gxxx(:,1,1,1),gxxx(:,2,1,1), ... 
                % because we have gxxx(beta1,gama1,gama2,gama3)
                i4 = 0;
            	for gama3=1:nx                  %for index controlling gxxx(:,:,:,here)
                    for gama2=1:nx              %for index controlling gxxx(:,:,here,:)
                        for gama1=1:nx          %for index controlling gxxx(:,here,:,:)
                            %for beta1=1:ny      %for index controlling gxxx(here,:,:,:)
                            %    i4 = i4+1;									
                            %    Q(m,i4) = fyp(i,beta1)*hx(gama3,alfa3)*hx(gama2,alfa2)*hx(gama1,alfa1);
                            %end
                            i4 = i4 + 1;
                            i1 = 1 + (i4-1)*ny;
                            i2 = i4*ny;
                            Q(m,i1:i2) = fyp(i,1:ny)*hx(gama3,alfa3)*hx(gama2,alfa2)*hx(gama1,alfa1);
                        end
                    end
                end

                % Q6 - adding to the constant a(:,1) and given loadings for hxxx                
                tmp = 0;
                for beta1=1:ny
                   tmp = tmp + (reshape(fypyp(i,beta1,:),1,ny)*gx_hx(:,alfa3)...
                               +reshape(fypy (i,beta1,:),1,ny)*gx(:,alfa3)...
                               +reshape(fypxp(i,beta1,:),1,nx)*hx(:,alfa3)...
                                       +fypx (i,beta1,alfa3))*gx(beta1,:)*squeeze(hxx(:,alfa1,alfa2))...
                         +fyp(i,beta1)*reshape(hxx(:,alfa1,alfa2),1,nx)*squeeze(gxx(beta1,:,:))*hx(:,alfa3);
                end
                a(m,1) = a(m,1) + tmp;
                % The loadings. Recall the first part of x is coefficient for gxxx, the second part is coefficients for hxxx.
                i1 = nxnxnxnyI+1+nx*(alfa1-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa3-1); 		
         		i2 = nxnxnxnyI+  nx *alfa1   +nxnxI*(alfa2-1)+nxnxnxI*(alfa3-1);
        	    Q(m,i1:i2) = fyp_gx(1,1:nx);    
                % The same as:  
                %i4 = nxnxnxnyI;
            	%for gama3=1:nx                  %for index controlling hxxx(:,:,:,here)
                %    for gama2=1:nx              %for index controlling hxxx(:,:,here,:)
                %        for gama1=1:nx          %for index controlling hxxx(:,here,:,:)
                %            for phi=1:nx        %for index controlling hxxx(here,:,:,:)
                %                i4 = i4+1;				
                %                if alfa1==gama1 && alfa2==gama2 && alfa3==gama3
                %                    Q(m,i4) = fyp_gx(1,phi);
                %                end
                %            end
                %       end
                %    end
                %end
                
                
                % Q7 - adding to the constant a(:,1)
                %tmp = 0;
                %for beta1=1:ny
                %    for beta2=1:ny
                %        tmp = tmp + (reshape(fyypyp(i,beta1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                %                    +reshape(fyypy (i,beta1,beta2,:),1,ny)*gx(:,alfa3)...
                %                    +reshape(fyypxp(i,beta1,beta2,:),1,nx)*hx(:,alfa3)...
                %                            +fyypx (i,beta1,beta2,alfa3))*gx_hx(beta2,alfa2)*gx(beta1,alfa1)...
                %              +fyyp(i,beta1,beta2)*(hx(:,alfa2)'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)*gx(beta1,alfa1)...
                %                    +gx(beta2,:)*squeeze(hxx(:,alfa2,alfa3))*gx(beta1,alfa1)...
                %                    +gx_hx(beta2,alfa2)*gxx(beta1,alfa1,alfa3));
                %    end
                %end
                tmp = 0;
                for beta2=1:ny
                    tmp = tmp + gx(:,alfa1)'*(reshape(fyypyp(i,:,beta2,:),ny,ny)*gx_hx(:,alfa3)...
                                +reshape(fyypy (i,:,beta2,:),ny,ny)*gx(:,alfa3)...
                                +reshape(fyypxp(i,:,beta2,:),ny,nx)*hx(:,alfa3)...
                                +reshape(fyypx (i,:,beta2,alfa3),ny,1))*gx_hx(beta2,alfa2)...
                          +reshape(fyyp(i,:,beta2),1,ny)*gx(:,alfa1)...
                          *hx(:,alfa2)'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)...
                          +reshape(fyyp(i,:,beta2),1,ny)*gx(:,alfa1)...
                          *gx(beta2,:)*squeeze(hxx(:,alfa2,alfa3))...
                          +reshape(fyyp(i,:,beta2),1,ny)*squeeze(gxx(:,alfa1,alfa3))...
                          *gx_hx(beta2,alfa2);
                end
                a(m,1) = a(m,1) + tmp;
                
                
                % Q8 - adding to the constant a(:,1)
                %tmp = 0;
                %for beta1=1:ny
                %    for beta2=1:ny
                %        tmp = tmp + (reshape(fyyyp(i,beta1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                %                    +reshape(fyyy (i,beta1,beta2,:),1,ny)*gx(:,alfa3)...
                %                    +reshape(fyyxp(i,beta1,beta2,:),1,nx)*hx(:,alfa3)...
                %                            +fyyx (i,beta1,beta2,alfa3))*gx(beta2,alfa2)*gx(beta1,alfa1)...
                %              +fyy(i,beta1,beta2)*(gxx(beta2,alfa2,alfa3)*gx(beta1,alfa1)...
                %                    +gx(beta2,alfa2)*gxx(beta1,alfa1,alfa3));
                %    end
                %end
                tmp = 0;
                for beta2=1:ny
                     tmp = tmp + gx(:,alfa1)'*(reshape(fyyyp(i,:,beta2,:),ny,ny)*gx_hx(:,alfa3)...
                                 +reshape(fyyy (i,:,beta2,:),ny,ny)*gx(:,alfa3)...
                                 +reshape(fyyxp(i,:,beta2,:),ny,nx)*hx(:,alfa3)...
                                 +reshape(fyyx (i,:,beta2,alfa3),ny,1))*gx(beta2,alfa2)...
                           +gx(:,alfa1)'*reshape(fyy(i,:,beta2),ny,1)*gxx(beta2,alfa2,alfa3)...
                           +reshape(fyy(i,:,beta2),1,ny)*squeeze(gxx(:,alfa1,alfa3))*gx(beta2,alfa2);
                end
                a(m,1) = a(m,1) + tmp;


                % Q9 - adding to the constant a(:,1)
                %tmp = 0;
                %for beta1=1:ny
                %    for gama2=1:nx
                %        tmp = tmp + (reshape(fyxpyp(i,beta1,gama2,:),1,ny)*gx_hx(:,alfa3)...
                %                    +reshape(fyxpy (i,beta1,gama2,:),1,ny)*gx(:,alfa3)...
                %                    +reshape(fyxpxp(i,beta1,gama2,:),1,nx)*hx(:,alfa3)...
                %                            +fyxpx (i,beta1,gama2,alfa3)*hx(gama2,alfa2)*gx(beta1,alfa1)...
                %              +fyxp(i,beta1,gama2)*(hxx(gama2,alfa2,alfa3)*gx(beta1,alfa1)...
                %                    +hx(gama2,alfa2)*gxx(beta1,alfa1,alfa3));
                %    end
                %end
                tmp = 0;
                for beta1=1:ny
                    tmp = tmp + hx(:,alfa2)'*(reshape(fyxpyp(i,beta1,:,:),nx,ny)*gx_hx(:,alfa3)...
                                +reshape(fyxpy (i,beta1,:,:),nx,ny)*gx(:,alfa3)...
                                +reshape(fyxpxp(i,beta1,:,:),nx,nx)*hx(:,alfa3)...
                                +reshape(fyxpx (i,beta1,:,alfa3),nx,1))*gx(beta1,alfa1)...
                          +reshape(fyxp(i,beta1,:),1,nx)*squeeze(hxx(:,alfa2,alfa3))*gx(beta1,alfa1)...
                          +reshape(fyxp(i,beta1,:),1,nx)*hx(:,alfa2)*gxx(beta1,alfa1,alfa3);
                end
                a(m,1) = a(m,1) + tmp;
                

                % Q10 - adding to the constant a(:,1)
                %tmp = 0;
                %for beta1=1:ny
                %    tmp = tmp + (reshape(fyxyp(i,beta1,alfa2,:),1,ny)*gx_hx(:,alfa3)...
                %                +reshape(fyxy (i,beta1,alfa2,:),1,ny)*gx(:,alfa3)...
                %                +reshape(fyxxp(i,beta1,alfa2,:),1,nx)*hx(:,alfa3)...
                %                        +fyxx (i,beta1,alfa2,alfa3))*gx(beta1,alfa1)...
                %          +fyx(i,beta1,alfa2)*gxx(beta1,alfa1,alfa3);
                %end
                %a(m,1) = a(m,1) + tmp;
                vec_ny1(1:ny,1)= reshape(fyxyp(i,:,alfa2,:),ny,ny)*gx_hx(:,alfa3)...
                                +reshape(fyxy (i,:,alfa2,:),ny,ny)*gx(:,alfa3)...
                                +reshape(fyxxp(i,:,alfa2,:),ny,nx)*hx(:,alfa3)...
                                +reshape(fyxx (i,:,alfa2,alfa3),ny,1);
                a(m,1) = a(m,1) + vec_ny1(1:ny,1)'*gx(1:ny,alfa1)+...
                                + reshape(fyx(i,:,alfa2),1,ny)*reshape(gxx(:,alfa1,alfa3),ny,1);
                
                
                % Q11 - adding to the constant a(:,1) and adding to the loadings for gxxx
                %tmp = 0;
                %for beta1=1:ny
                %   tmp = tmp + (reshape(fyyp(i,beta1,:),1,ny)*gx_hx(:,alfa3)...
                %               +reshape(fyy (i,beta1,:),1,ny)*gx(:,alfa3)...
                %               +reshape(fyxp(i,beta1,:),1,nx)*hx(:,alfa3)...                               
                %                       +fyx (i,beta1,alfa3))*gxx(beta1,alfa1,alfa2);
                %end
                %a(m,1) = a(m,1) + tmp;
                vec_ny(1:ny,1) = reshape(fyyp(i,:,:),ny,ny)*gx_hx(:,alfa3)...
                                +reshape(fyy (i,:,:),ny,ny)*gx(:,alfa3)...
                                +reshape(fyxp(i,:,:),ny,nx)*hx(:,alfa3)...                               
                                +reshape(fyx (i,:,alfa3),ny,1);
                a(m,1) = a(m,1) + vec_ny(1:ny,1)'*squeeze(gxx(:,alfa1,alfa2));
                % The loadings for gxxx
                i1				= 1+ny*(alfa1-1)+nxnyI*(alfa2-1)+nxnxnyI*(alfa3-1); 	
        		i2				=   ny* alfa1   +nxnyI*(alfa2-1)+nxnxnyI*(alfa3-1);	
        		Q(m,i1:i2)	    = Q(m,i1:i2)+fy(i,1:ny);                
                %i4 = 0;
            	%for gama3=1:nx                  %for index controlling gxxx(:,:,:,here)
                %    for gama2=1:nx              %for index controlling gxxx(:,:,here,:)
                %        for gama1=1:nx          %for index controlling gxxx(:,here,:,:)
                %           for beta1=1:ny       %for index controlling gxxx(here,:,:,:)
                %                i4 = i4+1;
                %                if alfa1==gama1 && alfa2==gama2 && alfa3==gama3
                %                    Q(m,i4) = Q(m,i4) + fy(i,beta1);
                %                end
                %            end
                %        end
                %    end
                %end                
            

                % Q12 - adding to the constant a(:,1)                
                %tmp = 0;
                %for gama1=1:nx
                %    for beta2=1:ny
                %        tmp = tmp + (reshape(fxpypyp(i,gama1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                %                    +reshape(fxpypy (i,gama1,beta2,:),1,ny)*gx(:,alfa3)...
                %                    +reshape(fxpypxp(i,gama1,beta2,:),1,nx)*hx(:,alfa3)...
                %                            +fxpypx (i,gama1,beta2,alfa3))*gx_hx(beta2,alfa2)*hx(gama1,alfa1)...
                %              +fxpyp(i,gama1,beta2)*(hx(:,alfa2)'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)*hx(gama1,alfa1)...
                %                    +gx(beta2,:)*squeeze(hxx(:,alfa2,alfa3))*hx(gama1,alfa1)...
                %                    +gx_hx(beta2,alfa2)*hxx(gama1,alfa1,alfa3));
                %    end
                %end
                tmp = 0;
                for beta2=1:ny
                     tmp = tmp + hx(:,alfa1)'*(reshape(fxpypyp(i,:,beta2,:),nx,ny)*gx_hx(:,alfa3)...
                                 +reshape(fxpypy (i,:,beta2,:),nx,ny)*gx(:,alfa3)...
                                 +reshape(fxpypxp(i,:,beta2,:),nx,nx)*hx(:,alfa3)...
                                 +reshape(fxpypx(i,:,beta2,alfa3),nx,1))*gx_hx(beta2,alfa2)...
                           +reshape(fxpyp(i,:,beta2),1,nx)*hx(:,alfa1)...
                           *hx(:,alfa2)'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)...
                           +reshape(fxpyp(i,:,beta2),1,nx)*hx(:,alfa1)...
                           *gx(beta2,:)*squeeze(hxx(:,alfa2,alfa3))...
                           +reshape(fxpyp(i,:,beta2),1,nx)*squeeze(hxx(:,alfa1,alfa3))*gx_hx(beta2,alfa2);
                end
                a(m,1) = a(m,1) + tmp;

                
                % Q13 - adding to the constant a(:,1)      
                %tmp = 0;
                %for gama1=1:nx
                %    for beta2=1:ny
                %        tmp = tmp + (reshape(fxpyyp(i,gama1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                %                    +reshape(fxpyy (i,gama1,beta2,:),1,ny)*gx(:,alfa3)...
                %                    +reshape(fxpyxp(i,gama1,beta2,:),1,nx)*hx(:,alfa3)...
                %                            +fxpyx (i,gama1,beta2,alfa3))*gx(beta2,alfa2)*hx(gama1,alfa1)...
                %              +fxpy(i,gama1,beta2)*(gxx(beta2,alfa2,alfa3)*hx(gama1,alfa1)...
                %                    +gx(beta2,alfa2)*hxx(gama1,alfa1,alfa3));
                %    end
                %end
                tmp = 0;
                for beta2=1:ny
                    tmp = tmp + hx(:,alfa1)'*(reshape(fxpyyp(i,:,beta2,:),nx,ny)*gx_hx(:,alfa3)...
                                +reshape(fxpyy (i,:,beta2,:),nx,ny)*gx(:,alfa3)...
                                +reshape(fxpyxp(i,:,beta2,:),nx,nx)*hx(:,alfa3)...
                                +reshape(fxpyx (i,:,beta2,alfa3),nx,1))*gx(beta2,alfa2)...
                          +reshape(fxpy(i,:,beta2),1,nx)*hx(:,alfa1)*gxx(beta2,alfa2,alfa3)...
                          +reshape(fxpy(i,:,beta2),1,nx)*squeeze(hxx(:,alfa1,alfa3))*gx(beta2,alfa2);
                end
                a(m,1) = a(m,1) + tmp;
                
                
                % Q14 - adding to the constant a(:,1)                      
                %tmp = 0;
                %for gama1=1:nx
                %    for gama2=1:nx
                %        tmp = tmp + (reshape(fxpxpyp(i,gama1,gama2,:),1,ny)*gx_hx(:,alfa3)...
                %                    +reshape(fxpxpy (i,gama1,gama2,:),1,ny)*gx(:,alfa3)...
                %                    +reshape(fxpxpxp(i,gama1,gama2,:),1,nx)*hx(:,alfa3)...
                %                            +fxpxpx (i,gama1,gama2,alfa3))*hx(gama2,alfa2)*hx(gama1,alfa1)...
                %              +fxpxp(i,gama1,gama2)*(hxx(gama2,alfa2,alfa3)*hx(gama1,alfa1)...
                %                    +hx(gama2,alfa2)*hxx(gama1,alfa1,alfa3));
                %    end
                %end
                tmp = 0;
                for gama2=1:nx
                    tmp = tmp + hx(:,alfa1)'*(reshape(fxpxpyp(i,:,gama2,:),nx,ny)*gx_hx(:,alfa3)...
                                +reshape(fxpxpy (i,:,gama2,:),nx,ny)*gx(:,alfa3)...
                                +reshape(fxpxpxp(i,:,gama2,:),nx,nx)*hx(:,alfa3)...
                                +reshape(fxpxpx (i,:,gama2,alfa3),nx,1))*hx(gama2,alfa2)...
                          +reshape(fxpxp(i,:,gama2),1,nx)*hx(:,alfa1)*hxx(gama2,alfa2,alfa3)...
                          +reshape(fxpxp(i,:,gama2),1,nx)*squeeze(hxx(:,alfa1,alfa3))*hx(gama2,alfa2);
                end
                a(m,1) = a(m,1) + tmp;

                
                % Q15 - adding to the constant a(:,1)                   
                %tmp = 0;
                %for gama1=1:nx
                %    tmp = tmp +  (reshape(fxpxyp(i,gama1,alfa2,:),1,ny)*gx_hx(:,alfa3)...
                %                 +reshape(fxpxy (i,gama1,alfa2,:),1,ny)*gx(:,alfa3)...
                %                 +reshape(fxpxxp(i,gama1,alfa2,:),1,nx)*hx(:,alfa3)...
                %                         +fxpxx (i,gama1,alfa2,alfa3))*hx(gama1,alfa1)...
                %              +fxpx(i,gama1,alfa2)*hxx(gama1,alfa1,alfa3);
                %end
                %a(m,1) = a(m,1) + tmp;
                vec_nx(1:nx,1) = reshape(fxpxyp(i,:,alfa2,:),nx,ny)*gx_hx(:,alfa3)...
                                +reshape(fxpxy (i,:,alfa2,:),nx,ny)*gx(:,alfa3)...
                                +reshape(fxpxxp(i,:,alfa2,:),nx,nx)*hx(:,alfa3)...
                                +reshape(fxpxx(i,:,alfa2,alfa3),nx,1);
                a(m,1) = a(m,1) + vec_nx(1:nx,1)'*hx(:,alfa1)...
                                + reshape(fxpx(i,:,alfa2),1,nx)*squeeze(hxx(:,alfa1,alfa3));
                
                
                % Q16 - adding to the constant a(:,1) and the loadings for hxxx
                %tmp = 0;
                %for gama1=1:nx
                %    tmp = tmp + (reshape(fxpyp(i,gama1,:),1,ny)*gx_hx(:,alfa3)...
                %                +reshape(fxpy (i,gama1,:),1,ny)*gx(:,alfa3)...
                %                +reshape(fxpxp(i,gama1,:),1,nx)*hx(:,alfa3)...
                %                        +fxpx (i,gama1,alfa3))*hxx(gama1,alfa1,alfa2);
                %end
                %a(m,1) = a(m,1) + tmp;
                vec_nx(1:nx,1) = reshape(fxpyp(i,:,:),nx,ny)*gx_hx(:,alfa3)...
                                +reshape(fxpy (i,:,:),nx,ny)*gx(:,alfa3)...
                                +reshape(fxpxp(i,:,:),nx,nx)*hx(:,alfa3)...
                                +reshape(fxpx (i,:,alfa3),nx,1);
                a(m,1) = a(m,1) +vec_nx(1:nx,1)'*squeeze(hxx(:,alfa1,alfa2));
                                
                % The loadings for hxxx
                i1 = nxnxnxnyI+1+nx*(alfa1-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa3-1); 		
         		i2 = nxnxnxnyI+  nx* alfa1   +nxnxI*(alfa2-1)+nxnxnxI*(alfa3-1);
                Q(m,i1:i2) = Q(m,i1:i2) + fxp(i,1:nx);
                %The same as
                %i4 = nxnxnxnyI;
            	%for gama3=1:nx                  %for index controlling hxxx(:,:,:,here)
                %    for gama2=1:nx              %for index controlling hxxx(:,:,here,:)
                %        for gama1=1:nx          %for index controlling hxxx(:,here,:,:)
                %           for phi=1:nx        %for index controlling hxxx(here,:,:,:)
                %                i4 = i4+1;				
                %                if alfa1==gama1 && alfa2==gama2 && alfa3==gama3
                %                    Q(m,i4) = Q(m,i4) + fxp(i,phi);
                %               end
                %           end
                %        end
                %    end
                %end
                
                
                % Q17 - adding to the constant a(:,1)
                tmp = 0;
                for beta2=1:ny
                   tmp = tmp + (reshape(fxypyp(i,alfa1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                               +reshape(fxypy (i,alfa1,beta2,:),1,ny)*gx(:,alfa3)...
                               +reshape(fxypxp(i,alfa1,beta2,:),1,nx)*hx(:,alfa3)...
                                       +fxypx (i,alfa1,beta2,alfa3))*gx_hx(beta2,alfa2)...
                         +fxyp(i,alfa1,beta2)*(hx(:,alfa2)'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)...
                               +gx(beta2,:)*squeeze(hxx(:,alfa2,alfa3)));
                end
                a(m,1) = a(m,1) + tmp;
                
                
                % Q18 - adding to the constant a(:,1)
                %tmp = 0;
                %for beta2=1:ny
                %   tmp = tmp + (reshape(fxyyp(i,alfa1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                %               +reshape(fxyy (i,alfa1,beta2,:),1,ny)*gx(:,alfa3)...
                %               +reshape(fxyxp(i,alfa1,beta2,:),1,nx)*hx(:,alfa3)...
                %                       +fxyx (i,alfa1,beta2,alfa3))*gx(beta2,alfa2)...
                %         +fxy(i,alfa1,beta2)*gxx(beta2,alfa2,alfa3);
                %end
                %a(m,1) = a(m,1) + tmp;
                vec_ny(1:ny,1) = reshape(fxyyp(i,alfa1,:,:),ny,ny)*gx_hx(:,alfa3)...
                                +reshape(fxyy (i,alfa1,:,:),ny,ny)*gx(:,alfa3)...
                                +reshape(fxyxp(i,alfa1,:,:),ny,nx)*hx(:,alfa3)...
                                +reshape(fxyx(i,alfa1,:,alfa3),ny,1);
                a(m,1) = a(m,1) +vec_ny(1:ny,1)'*gx(:,alfa2)...
                                +reshape(fxy(i,alfa1,:),1,ny)*squeeze(gxx(:,alfa2,alfa3));

                            
                % Q19 - adding to the constant a(:,1)
                %tmp = 0;
                %for gama2=1:nx
                %   tmp = tmp + (reshape(fxxpyp(i,alfa1,gama2,:),1,ny)*gx_hx(:,alfa3)... 
                %               +reshape(fxxpy (i,alfa1,gama2,:),1,ny)*gx(:,alfa3)...
                %               +reshape(fxxpxp(i,alfa1,gama2,:),1,nx)*hx(:,alfa3)...
                %                       +fxxpx (i,alfa1,gama2,alfa3))*hx(gama2,alfa2)...
                %         +fxxp(i,alfa1,gama2)*hxx(gama2,alfa2,alfa3);
                %end
                %a(m,1) = a(m,1) + tmp;
                vec_nx(1:nx,1) = reshape(fxxpyp(i,alfa1,:,:),nx,ny)*gx_hx(:,alfa3)... 
                                +reshape(fxxpy (i,alfa1,:,:),nx,ny)*gx(:,alfa3)...
                                +reshape(fxxpxp(i,alfa1,:,:),nx,nx)*hx(:,alfa3)...
                                +reshape(fxxpx (i,alfa1,:,alfa3),nx,1);
                a(m,1) = a(m,1) +vec_nx(1:nx,1)'*hx(:,alfa2)...
                                +reshape(fxxp(i,alfa1,:),1,nx)*squeeze(hxx(:,alfa2,alfa3));
                
                    
                % Q20 - adding to the constant a(:,1)
                a(m,1) = a(m,1) + reshape(fxxyp(i,alfa1,alfa2,:),1,ny)*gx_hx(:,alfa3)... 
                                 +reshape(fxxy (i,alfa1,alfa2,:),1,ny)*gx(:,alfa3)...
                                 +reshape(fxxxp(i,alfa1,alfa2,:),1,nx)*hx(:,alfa3)...
                                 +fxxx (i,alfa1,alfa2,alfa3);                  
            end
        end
    end    
end

% Filling data into Qr(colQ,colQ). Here we reduce the numbers of columns in Q from
% n*nx*nx*nx to a colQ
my = 0;
mx = 0;
num_gxxx = ny*(nx^2+nx*(nx-1)*(nx-2)/6);
for alfa3=1:nx
    for alfa2=alfa3:nx
        for alfa1=alfa2:nx
            % The loadings for gxxx
            for i=1:ny
                my = my+1;
                if     alfa1 == alfa2 && alfa2 == alfa3 
                    Qr(1:colQ,my) = Q(1:colQ,ny*(alfa1-1)+nxnyI*(alfa2-1)+nxnxnyI*(alfa3-1)+i);
                elseif alfa1 == alfa2 && alfa2 ~= alfa3                    %alfa1==alfa2~=alfa3
                    Qr(1:colQ,my) = Q(1:colQ,ny*(alfa1-1)+nxnyI*(alfa1-1)+nxnxnyI*(alfa3-1)+i)...
                                   +Q(1:colQ,ny*(alfa1-1)+nxnyI*(alfa3-1)+nxnxnyI*(alfa1-1)+i)...
                                   +Q(1:colQ,ny*(alfa3-1)+nxnyI*(alfa1-1)+nxnxnyI*(alfa1-1)+i);
                elseif alfa1 ~= alfa2 && alfa2 == alfa3                    %alfa1~=alfa2==alfa3
                    Qr(1:colQ,my) = Q(1:colQ,ny*(alfa1-1)+nxnyI*(alfa2-1)+nxnxnyI*(alfa2-1)+i)...
                                   +Q(1:colQ,ny*(alfa2-1)+nxnyI*(alfa1-1)+nxnxnyI*(alfa2-1)+i)...
                                   +Q(1:colQ,ny*(alfa2-1)+nxnyI*(alfa2-1)+nxnxnyI*(alfa1-1)+i);
                elseif alfa1 ~= alfa2 && alfa1 ~= alfa3 &&  alfa2 ~= alfa3 %alfa1~=alfa2~=alfa3
                    Qr(1:colQ,my) = Q(1:colQ,ny*(alfa1-1)+nxnyI*(alfa2-1)+nxnxnyI*(alfa3-1)+i)...
                                   +Q(1:colQ,ny*(alfa1-1)+nxnyI*(alfa3-1)+nxnxnyI*(alfa2-1)+i)...
                                   +Q(1:colQ,ny*(alfa3-1)+nxnyI*(alfa1-1)+nxnxnyI*(alfa2-1)+i)...
                                   +Q(1:colQ,ny*(alfa3-1)+nxnyI*(alfa2-1)+nxnxnyI*(alfa1-1)+i)...
                                   +Q(1:colQ,ny*(alfa2-1)+nxnyI*(alfa3-1)+nxnxnyI*(alfa1-1)+i)...
                                   +Q(1:colQ,ny*(alfa2-1)+nxnyI*(alfa1-1)+nxnxnyI*(alfa3-1)+i);
                end
            end
            % The loadings for hxxx
            for i=1:nx
                mx = mx+1;
                if     alfa1 == alfa2 && alfa2 == alfa3 
                    Qr(1:colQ,mx+num_gxxx) = Q(1:colQ,nxnxnxnyI+nx*(alfa1-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa3-1)+i);
                elseif alfa1 == alfa2 && alfa2 ~= alfa3 %alfa1==alfa2~=alfa3
                    Qr(1:colQ,mx+num_gxxx) = ...
                                    Q(1:colQ,nxnxnxnyI+nx*(alfa1-1)+nxnxI*(alfa1-1)+nxnxnxI*(alfa3-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa1-1)+nxnxI*(alfa3-1)+nxnxnxI*(alfa1-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa3-1)+nxnxI*(alfa1-1)+nxnxnxI*(alfa1-1)+i);
                elseif alfa1 ~= alfa2 && alfa2 == alfa3 %alfa1~=alfa2==alfa3
                    Qr(1:colQ,mx+num_gxxx) = ...
                                    Q(1:colQ,nxnxnxnyI+nx*(alfa1-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa2-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa2-1)+nxnxI*(alfa1-1)+nxnxnxI*(alfa2-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa2-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa1-1)+i);                    
                               
                elseif alfa1 ~= alfa2 && alfa1 ~= alfa3 &&  alfa2 ~= alfa3 %alfa1~=alfa2~=alfa3
                    Qr(1:colQ,mx+num_gxxx) = ...
                                    Q(1:colQ,nxnxnxnyI+nx*(alfa1-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa3-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa1-1)+nxnxI*(alfa3-1)+nxnxnxI*(alfa2-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa3-1)+nxnxI*(alfa1-1)+nxnxnxI*(alfa2-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa3-1)+nxnxI*(alfa2-1)+nxnxnxI*(alfa1-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa2-1)+nxnxI*(alfa3-1)+nxnxnxI*(alfa1-1)+i)...
                                   +Q(1:colQ,nxnxnxnyI+nx*(alfa2-1)+nxnxI*(alfa1-1)+nxnxnxI*(alfa3-1)+i);
                end
            end
        end
    end
end



% we solve the system Qr(colQ,colQ)*x(colQ,1) + a(colQ,1) = 0
x = Qr\(-a(:,1)); %same as: inv(Qr)*(-a(:,1))

% Representing the solution for gxxx...
m = 0;
for alfa3=1:nx
    for alfa2=alfa3:nx
        for alfa1=alfa2:nx
            for i=1:ny
                m = m+1;
                gxxx(i,alfa1,alfa2,alfa3) = x(m);
            end
            % Using symmetry for alfa1 and alfa2
            if alfa1 == alfa2 && alfa2 ~= alfa3 %alfa1==alfa2~=alfa3
                %gxxx(:,alfa1,alfa1,alfa3)= gxxx(:,alfa1,alfa2,alfa3);
                gxxx(:,alfa1,alfa3,alfa1) = gxxx(:,alfa1,alfa2,alfa3);
                gxxx(:,alfa3,alfa1,alfa1) = gxxx(:,alfa1,alfa2,alfa3);                
            end
            % Using symmetry for alfa2 and alfa3            
            if alfa1 ~= alfa2 && alfa2 == alfa3  %alfa1~=alfa2==alfa3 
                %gxxx(:,alfa1,alfa2,alfa2)= gxxx(:,alfa1,alfa2,alfa3);                
                gxxx(:,alfa2,alfa1,alfa2) = gxxx(:,alfa1,alfa2,alfa3);                
                gxxx(:,alfa2,alfa2,alfa1) = gxxx(:,alfa1,alfa2,alfa3);                                
            end            
            % Using symmetry for alfa1,alfa2, and alfa3            
            if alfa1 ~= alfa2 && alfa1 ~= alfa3 &&  alfa2 ~= alfa3 %alfa1~=alfa2~=alfa3
                %gxxx(:,alfa1,alfa2,alfa3) = gxxx(:,alfa1,alfa2,alfa3);
                gxxx(:,alfa1,alfa3,alfa2) = gxxx(:,alfa1,alfa2,alfa3);
                gxxx(:,alfa3,alfa1,alfa2) = gxxx(:,alfa1,alfa2,alfa3);
                gxxx(:,alfa3,alfa2,alfa1) = gxxx(:,alfa1,alfa2,alfa3);
                gxxx(:,alfa2,alfa3,alfa1) = gxxx(:,alfa1,alfa2,alfa3); 
                gxxx(:,alfa2,alfa1,alfa3) = gxxx(:,alfa1,alfa2,alfa3);                 
            end
        end
    end
end

% ...and for hxxx
for alfa3=1:nx    
    for alfa2=alfa3:nx
        for alfa1=alfa2:nx
            for i=1:nx
                m = m+1;
                hxxx(i,alfa1,alfa2,alfa3) = x(m);
            end
           % Using symmetry for alfa1 and alfa2
            if alfa1 == alfa2 && alfa2 ~= alfa3 %alfa1==alfa2~=alfa3
                %hxxx(:,alfa1,alfa1,alfa3)= hxxx(:,alfa1,alfa2,alfa3);
                hxxx(:,alfa1,alfa3,alfa1) = hxxx(:,alfa1,alfa2,alfa3);
                hxxx(:,alfa3,alfa1,alfa1) = hxxx(:,alfa1,alfa2,alfa3);                
            end
            % Using symmetry for alfa2 and alfa3            
            if alfa1 ~= alfa2 && alfa2 == alfa3  %alfa1~=alfa2==alfa3 
                %hxxx(:,alfa1,alfa2,alfa2)= hxxx(:,alfa1,alfa2,alfa3);                
                hxxx(:,alfa2,alfa1,alfa2) = hxxx(:,alfa1,alfa2,alfa3);                
                hxxx(:,alfa2,alfa2,alfa1) = hxxx(:,alfa1,alfa2,alfa3);                                
            end            
            % Using symmetry for alfa1,alfa2, and alfa3            
            if alfa1 ~= alfa2 && alfa1 ~= alfa3 &&  alfa2 ~= alfa3 %alfa1~=alfa2~=alfa3
                %hxxx(:,alfa1,alfa2,alfa3) = hxxx(:,alfa1,alfa2,alfa3);
                hxxx(:,alfa1,alfa3,alfa2) = hxxx(:,alfa1,alfa2,alfa3);
                hxxx(:,alfa3,alfa1,alfa2) = hxxx(:,alfa1,alfa2,alfa3);
                hxxx(:,alfa3,alfa2,alfa1) = hxxx(:,alfa1,alfa2,alfa3);
                hxxx(:,alfa2,alfa3,alfa1) = hxxx(:,alfa1,alfa2,alfa3); 
                hxxx(:,alfa2,alfa1,alfa3) = hxxx(:,alfa1,alfa2,alfa3);                 
            end
        end
    end
end


% The system we solve to get gssx and hssx: Qss(nnxI,nnxI)*xss(nnxI,1) + ass(nnxI,1) = 0
% The index m counts the number of unknowns, i.e. the elements in x
gxEta(1:ny,1:ne) = gx*eta;
I_ne           = eye(ne);

m = 0;
for i=1:n
   for alfa3=1:nx
       m = m + 1;
       %1) Adding to ass
       tmp = 0;
       for beta1=1:ny
           for beta2=1:ny
               tmp = tmp + (reshape(fypypyp(i,beta1,beta2,:),1,ny)*gx_hx(:,alfa3)...
                           +reshape(fypypy (i,beta1,beta2,:),1,ny)*gx(:,alfa3)...
                           +reshape(fypypxp(i,beta1,beta2,:),1,nx)*hx(:,alfa3)...
                           +fypypx(i,beta1,beta2,alfa3))*gxEta(beta1,:)*I_ne*gxEta(beta2,:)'...
                     +fypyp(i,beta1,beta2)*(gxEta(beta1,:)*I_ne*eta'*squeeze(gxx(beta2,:,:))*hx(:,alfa3)...
                           + gxEta(beta2,:)*I_ne*eta'*squeeze(gxx(beta1,:,:))*hx(:,alfa3));
           end
       end
       ass(m,1) = tmp;
       

       %4) Adding to ass       
       tmp = 0;
       for beta1=1:ny
           for gama2=1:nx
               MM_gama1(1:nx,1) = eta(1:nx,:)*I_ne*eta(gama2,:)';
               tmp = tmp+2*(reshape(fypxpyp(i,beta1,gama2,:),1,ny)*gx_hx(:,alfa3)...
                           +reshape(fypxpy (i,beta1,gama2,:),1,ny)*gx(:,alfa3)...
                           +reshape(fypxpxp(i,beta1,gama2,:),1,nx)*hx(:,alfa3)...
                           +fypxpx(i,beta1,gama2,alfa3))*gxEta(beta1,:)*I_ne*eta(gama2,:)'...
                     +2*fypxp(i,beta1,gama2)*MM_gama1'*squeeze(gxx(beta1,:,:))*hx(:,alfa3);
           end
       end
       ass(m,1) = ass(m,1) + tmp;

       
       %5) Adding to ass 
       tmp  = 0;
       for beta1=1:ny
           tmp1 = 0;
           tmp2 = 0;           
           for gama1=1:nx
               for gama2=1:nx
                  tmp1 = tmp1 + gxx(beta1,gama1,gama2)*eta(gama2,:)*I_ne*eta(gama1,:)';
                  tmp2 = tmp2 + reshape(gxxx(beta1,gama1,gama2,:),1,nx)*hx(:,alfa3)...
                                *eta(gama2,:)*I_ne*eta(gama1,:)';
               end
           end
           tmp = tmp + (reshape(fypyp(i,beta1,:),1,ny)*gx_hx(:,alfa3)...
                       +reshape(fypy (i,beta1,:),1,ny)*gx(:,alfa3)...
                       +reshape(fypxp(i,beta1,:),1,nx)*hx(:,alfa3)...
                       +fypx(i,beta1,alfa3))*(tmp1+gx(beta1,:)*hss+gss(beta1,1))...
                 +fyp(i,beta1)*(tmp2 + hss'*squeeze(gxx(beta1,:,:))*hx(:,alfa3));
       end
       ass(m,1) = ass(m,1) + tmp;

       
       %7) Adding to the loadings for xss
       % The loadings
       i1 = nxnyI+1+nx*(alfa3-1); 		
       i2 = nxnyI+nx*alfa3;
       Qss(m,i1:i2) = fyp(i,:)*gx(:,:);
       
       
       %13) Adding to the loadings for xss       
       % The loadings
       i4 = 0;
       for gama3=1:nx
           for beta1=1:ny
               i4 = i4 + 1;
               Qss(m,i4) = fyp(i,beta1)*hx(gama3,alfa3);
           end
       end       

       
       %18) Adding to ass and the loadings for xss           
       %tmp = 0;
       %for beta1=1:ny
       %    tmp = tmp + (reshape(fyyp(i,beta1,:),1,ny)*gx_hx(:,alfa3)...
       %                +reshape(fyy (i,beta1,:),1,ny)*gx(:,alfa3)...
       %                +reshape(fyxp(i,beta1,:),1,nx)*hx(:,alfa3)...
       %                +fyx(i,beta1,alfa3))*gss(beta1,1);
       %end
       vec_ny(1:ny,1) = reshape(fyyp(i,:,:),ny,ny)*gx_hx(:,alfa3)...
                       +reshape(fyy (i,:,:),ny,ny)*gx(:,alfa3)...
                       +reshape(fyxp(i,:,:),ny,nx)*hx(:,alfa3)...
                       +reshape(fyx (i,:,alfa3),ny,1);
       ass(m,1) = ass(m,1) + vec_ny(1:ny,1)'*gss(:,1);
       % The loadings
       i1 = 1+ny*(alfa3-1); 		
       i2 =   ny*alfa3;
       Qss(m,i1:i2) = Qss(m,i1:i2) + fy(i,1:ny);
     
       %22) Adding to ass        
       tmp = 0;
       for gama1=1:nx
           for gama2=1:nx
               tmp = tmp + (reshape(fxpxpyp(i,gama1,gama2,:),1,ny)*gx_hx(:,alfa3)...
                           +reshape(fxpxpy (i,gama1,gama2,:),1,ny)*gx(:,alfa3)...
                           +reshape(fxpxpxp(i,gama1,gama2,:),1,nx)*hx(:,alfa3)...
                           +fxpxpx(i,gama1,gama2,alfa3))*eta(gama1,:)*I_ne*eta(gama2,:)';
           end
       end
       ass(m,1) = ass(m,1) + tmp;          
       
   
       %23) Adding to ass and the loadings for xss                  
       %tmp = 0;
       %for gama1=1:nx
       %    tmp = tmp + (reshape(fxpyp(i,gama1,:),1,ny)*gx_hx(:,alfa3)...
       %                +reshape(fxpy (i,gama1,:),1,ny)*gx(:,alfa3)...
       %                +reshape(fxpxp(i,gama1,:),1,nx)*hx(:,alfa3)...
       %                +fxpx(i,gama1,alfa3))*hss(gama1,1);
       %end
       %ass(m,1) = ass(m,1) + tmp;   
       vec_nx(1:nx,1) = reshape(fxpyp(i,:,:),nx,ny)*gx_hx(:,alfa3)...
                       +reshape(fxpy (i,:,:),nx,ny)*gx(:,alfa3)...
                       +reshape(fxpxp(i,:,:),nx,nx)*hx(:,alfa3)...
                       +reshape(fxpx (i,:,alfa3),nx,1);
       ass(m,1) = ass(m,1) + vec_nx(1:nx,1)'*hss(:,1);
       % The loadings
       i1 = nxnyI+1+nx*(alfa3-1); 		
       i2 = nxnyI+  nx*alfa3;
       Qss(m,i1:i2) = Qss(m,i1:i2) + fxp(i,1:nx);       
   end
end

% we solve the system Qss(nnxI,nnxI)*xss(nnxI,1) + ass(nnxI,1) = 0
xss = Qss\(-ass(:,1)); %same as: inv(Qss)*(-ass(:,1)); 

% Representing the solution
for alfa3=1:nx
    i1 = 1+ny*(alfa3-1); 		
    i2 =   ny*alfa3;
    gssx(1:ny,alfa3) = xss(i1:i2);
    
    i1 = nxnyI+1+nx*(alfa3-1); 		
    i2 = nxnyI+  nx*alfa3;
    hssx(1:nx,alfa3) = xss(i1:i2);    
end

% If all third moments are zero, then
if any(vectorMom3) == 0
    gsss = zeros(ny,1);
    hsss = zeros(nx,1);
    return
end

% The system we solve to get gsss and hsss: Qsss(n,n)*xsss(n,1) + asss(n,1) = 0
% The index m counts the number of unknowns, i.e. the elements in xsss
Qsss = zeros(n,n);
asss = zeros(n,1);
% The value of asss
for i=1:n
    % Computing terms 1), 2), 3), 4) and 5)
    tmp = 0;
    for beta1=1:ny
        for phi1=1:ne
            tmp = tmp + (gx*eta(:,phi1))'*reshape(fypypyp(i,beta1,:,:),ny,ny)...                                 %1)
                             *gx*eta(:,phi1)*gx(beta1,:)*eta(:,phi1)*vectorMom3(phi1)...  
                      + 3*(gx*eta(:,phi1))'*reshape(fypypxp(i,beta1,:,:),ny,nx)...                               %2)
                             *eta(:,phi1)*gx(beta1,:)*eta(:,phi1)*vectorMom3(phi1)....
                      + 3*(gx(:,:)*eta(:,phi1)*vectorMom3(phi1))'*reshape(fypyp(i,:,beta1),ny,1)*eta(:,phi1)'... %3)
                             *squeeze(gxx(beta1,:,:))*eta(:,phi1)...
                      + 3*(eta(:,phi1))'*reshape(fypxpxp(i,beta1,:,:),nx,nx)...                                  %4)
                             *eta(:,phi1)*gx(beta1,:)*eta(:,phi1)*vectorMom3(phi1)...
                      + 3*reshape(fypxp(i,beta1,:),1,nx)*eta(:,phi1)...                                          %5)
                             *(eta(:,phi1)*vectorMom3(phi1))'*squeeze(gxx(beta1,:,:))*eta(:,phi1);
        end
    end
    asss(i,1) = tmp;

    % Computing term 6)
    tmp = 0;
    for beta1=1:ny
        for gama1=1:nx
             for phi1=1:ne
                    tmp = tmp + fyp(i,beta1)*eta(:,phi1)'*reshape(gxxx(beta1,gama1,:,:),nx,nx)*eta(:,phi1)...     %6)
                                *eta(gama1,phi1)*vectorMom3(phi1);
             end
        end
    end
    asss(i,1) = asss(i,1) + tmp;  
 
    % Computing term 7)
    tmp = 0;
    for gama1=1:nx
        for phi1=1:ne
             tmp = tmp + eta(:,phi1)'*reshape(fxpxpxp(i,gama1,:,:),nx,nx)*eta(:,phi1)...                     %7)
                        *eta(gama1,phi1)*vectorMom3(phi1);
        end
    end
    asss(i,1) = asss(i,1) + tmp;     
       
    % The loading for the matrix
    Qsss(i,1:ny)   = reshape(fyp(i,:),1,ny) + reshape(fy(i,:),1,ny);                                              
    Qsss(i,ny+1:n) = reshape(fyp(i,:),1,ny)*gx + reshape(fxp(i,:),1,nx);                                          
end

% we solve the system Qsss*xsss + asss = 0
xsss = Qsss\(-asss(:,1)); %same as: inv(Qsss)*(-asss(:,1)); 
gsss(1:ny,1) = xsss(1:ny);
hsss(1:nx,1) = xsss(ny+1:n);
