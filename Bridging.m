clc;

% Parameters:
% n = Dimension
% N = Length of the Frame
% L = The Erasure Set
% W = The Bridge Set

n = 2000;
N = 3000;
L = [1:100];
W = [2001:2100];

% The columns of F are a Gaussian randomly generated frame.
% The columns of G are the standard dual to F.

F = (1/sqrt(n))*randn(n,N);
S = F * F';
G = S \ F;

% % To create a Parseval frame comment out lines 17-19
% % and uncomment lines 24-27 below.
% 
% F = randn(N,n);
% [F,~] = qr(F,0);
% F = F';
% G = F;

% f is a random vector that we will try to recover 
% from frame coefficient erasures.

f = rand(n,1);
f = f ./ norm(f,2);

% FC are the frame coefficients of f.

FC = G' * f;

% We erase the frame coefficients indexed by
% L, 

FC(L) = zeros(size(L'));

% We compute f_R.

f_R = F * FC;

% We reconstruct the erased frame coefficients.

FRCL = G(:,L)' * f_R;
FRCB = G(:,W)' * f_R;

C = (F(:,L)'*G(:,W))\(F(:,L)'*G(:,L));

FC(L) = C' * (FC(W) - FRCB) + FRCL;

% We reconstruct the signal.

g = f_R + F(:,L) * FC(L);

% We compute the \ell^2 norm of the reconstruction
% error.

norm(f-g,2)
