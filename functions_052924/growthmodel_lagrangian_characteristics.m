function [dydt]=growthmodel_lagrangian_characteristics(t,y,J_cubic,tau)

% variables
N = numel(y)/2;
rho = y(1:N);
rho(rho<=0) = 0;

if any(rho<0)
    error('negative')
end

x = y(N+1:2*N);
x(x<=0) = 0;

% parameters
% tau = 96; % decay time for ligand
Factor = 0.3595; % dv/dx vs ERKi.^alpha
% alpha = 2.4097; % G(u) vs E(u)
alpha = 2.3981; % G(u) vs E(u)


% source
J = J_cubic(x(end)-x,repmat(t,size(x)));
J(J<=0) = 0;

% lambda_source=.2; % 1/3./(hx), renormalize to number of dx
% tau_source = 1*24;
% J = exp(-t/tau_source)*exp(-(x(end)-x)/lambda_source); % initialize source




% equations

drhodt = J - 1/tau.*rho - Factor.*rho.^alpha.*rho;
dxdt = Factor.*cumtrapz(x,rho.^alpha);

% drhodt(rho<=0) = 0;
% dxdt(x<0) = 0;


% put back to derivative
dydt = [drhodt;dxdt];




end

