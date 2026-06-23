function [dydt]=growthmodel_lagrangian(t,y,time_dependent)

% variables
phi = y(1);
xi = y(2);

% parameters
a1 = 0.05088; % dLdt/L vs total fractionGEM
% a2 = 7.06473694; % G(u) vs E(u)
a2 = 7.19766154; % G(u) vs E(u)
% a2 = 4.43;

A0 = 0.35998; % A(phi)

% f(u) 
if time_dependent % from linear model Erk ~ 1 + phi_ + u:phi_ + phi_:dpa + u:phi_:dpa, averaged data
%     af = (0.16677+t./24.*0.034583)./A0; % f(u)
%     bf = (0.30788+t./24.*(-0.043484))./A0; % f(u)
    af = (0.17355+t./24.*0.032759)./A0; % f(u)
    bf = (0.32718+t./24.*(-0.023222))./A0; % f(u)
else
    % af = 0.68059; % f(u)
    % bf = 0.66059; % f(u)
    af = 1.2; % f(u)
    bf = 0.4; % f(u)
end


% alpha = 2.4097; % G(u) vs E(u)
alpha = 2.39807; % G(u) vs E(u)

% alpha = 1; a2  = 1.2367;
% alpha = 2; a2  = 4.2968;
% alpha = 3; a2  = 15.8853;
% alpha = 4; a2  = 60.0677;


% alpha = 2.2185; % timeaveraged total%GEM vs <E>
Cf = 1/af*1/(alpha+1)*((af+bf)^(alpha+1)-bf^(alpha+1)); % integral of f(u)^alpha
% a2 = 5.8498./Cf; % timeaveraged total%GEM vs <E>

gamma = a1*a2.*A0^alpha.*Cf;


% equations
dphidt = gamma.*(1-phi).^alpha.*phi;
% dphidt = gamma.*(1-phi).^alpha.*phi+0.00035; % add arbitrary indefinite growth
dxidt = a1*a2*A0^alpha.*(1-phi).^alpha.*phi.*...
    1/af*1/(alpha+1).*((af.*(xi./phi)+bf).^(alpha+1)-bf^(alpha+1));



% put back to derivative
dydt = [dphidt;dxidt];

end

