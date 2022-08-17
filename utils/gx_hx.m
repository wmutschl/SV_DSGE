function [gx,hx,ERROR] = gx_hx(fy,fx,fyp,fxp);
%This program computes the matrices gx and hx that define the first-order approximation 
%of the DSGE model. That is, if 
%E_t[f(yp,y,xp,x)=0, then the solution is of the form
%xp = h(x,sigma) + sigma * eta * ep
%y = g(x,sigma).
%The first-order approximations to the functions g and h around the point (x,sigma)=(xbar,0), where xbar=h(xbar,0), are:
%h(x,sigma) = xbar + hx (x-xbar) 
%and
%g(x,sigma) = ybar + gx * (x-xbar),
%where ybar=g(xbar,0). 
%Inputs: fy fyp fx fxp
%Outputs: gx hx
%Calls solab.m (by Paul Klein)
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001
A = [-fxp -fyp];
B = [fx fy];

[gx,hx,ERROR]=solab(A,B,size(fx,2));

% Function: solab
%
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number nk
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% nk is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
% Calls: qzdiv
%
function [f,p,ERROR] = solab(a,b,nk);
ERROR = 0;
[s,t,q,z] = qz(a,b);            % upper triangular factorization of the matrix pencil b-za
[s,t,q,z] = qzdiv(1,s,t,q,z);   % reordering of generalized eigenvalues in ascending order
z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);

if any(any(isnan(z11),1),2) ~= 0 | any(any(isnan(z11),1),2) ~= 0
    f = NaN;
    p = NaN;
    ERROR = 1;
    return;
end

if rank(z11)<nk;
   %warning('Invertibility condition violated');
   f = NaN;
   p = NaN;
   ERROR = 1;
   return;
end
z11i = z11\eye(nk);
s11 = s(1:nk,1:nk);
t11 = t(1:nk,1:nk);

if abs(t(nk,nk))>abs(s(nk,nk)) | abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
   %warning('Wrong number of stable eigenvalues.') ;
   f = NaN;
   p = NaN;
   ERROR = 2;
   return;
end
dyn = s11\t11;
f = real(z21*z11i);
p = real(z11*dyn*z11i);


function [A,B,Q,Z] = qzdiv(stake,A,B,Q,Z)
%
% Written by Chris Sims
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right 
% corner, while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.
%
[n jnk] = size(A);
root = abs([diag(A) diag(B)]);
root(:,1) = root(:,1)-(root(:,1)<1.e-13).*(root(:,1)+root(:,2));
root(:,2) = root(:,2)./root(:,1);
for i = n:-1:1
   m=0;
   for j=i:-1:1
      if (root(j,2) > stake | root(j,2) < -.1) 
         m=j;
         break
      end
   end
   if (m==0) 
      return 
   end
   for k=m:1:i-1
      [A B Q Z] = qzswitch(k,A,B,Q,Z);
      tmp = root(k,2);
      root(k,2) = root(k+1,2);
      root(k+1,2) = tmp;
   end
end         


function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
%function [A,B,Q,Z] = qzswitch(i,A,B,Q,Z)
% Written by Chris Sims
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, interchanges
% diagonal elements i and i+1 of both A and B, while maintaining 
% Q'AZ' and Q'BZ' unchanged.  Does nothing if ratios of diagonal elements
% in A and B at i and i+1 are the same.  Aborts if diagonal elements of
% both A and B are zero at either position.
%
a = A(i,i); d = B(i,i); b = A(i,i+1); e = B(i,i+1);
c = A(i+1,i+1); f = B(i+1,i+1); 
wz = [c*e-f*b, (c*d-f*a)'];
xy = [(b*d-e*a)', (c*d-f*a)'];
n = sqrt(wz*wz');
m = sqrt(xy*xy');
if n == 0
   return
else
   wz = n\wz;
   xy = m\xy;
   wz = [wz; -wz(2)', wz(1)'];
   xy = [xy;-xy(2)', xy(1)'];
   A(i:i+1,:) = xy*A(i:i+1,:);
   B(i:i+1,:) = xy*B(i:i+1,:);
   A(:,i:i+1) = A(:,i:i+1)*wz;
   B(:,i:i+1) = B(:,i:i+1)*wz;
   Z(:,i:i+1) = Z(:,i:i+1)*wz;
   Q(i:i+1,:) = xy*Q(i:i+1,:);
end
