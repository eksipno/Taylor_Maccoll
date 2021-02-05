clear
close all
clc

gamma = 1.4;    % specific heat ratio
M = [1.1,1.25,1.5,2,4,8];  % mach number
in = numel(M);
jn = 200;
tmp = linspace(0,1,jn).^2;
beta = zeros(in,jn);
tc = beta;
for i=1:in
    betaMin = asin(1/M(i));    % minimum shock angle
    beta(i,:) = tmp*(pi/2-betaMin)+betaMin; % betaMin<=beta<=pi/2,non uniform grid
    for j=1:jn
        tc(i,j) = coneAngle(gamma,M(i),beta(i,j));  % cone angle:theta_c
        fprintf('%d/%d\n',i*j,in*jn);
    end   
    [~,I] = max(tc(i,:));   % find index
    plot(rad2deg(tc(i,1:I)),rad2deg(beta(i,1:I)),'.-','DisplayName',sprintf('M=%3.2f',M(i)));
    hold on
end
grid on
xlabel('cone angle \theta_c');
ylabel('shock angle \beta');
legend('Location','southeast')


%% functions
function tc = coneAngle(gamma,M,beta)
delta = obliqueShock(gamma,M,beta);
% Mn1 = M*sin(beta);
% Mn2 = sqrt((Mn1^2+(2/(gamma-1)))/(2*gamma/(gamma-1)*Mn1^2-1));
% M2 = Mn2/(sin(beta-delta));
% V = (1+2/((gamma-1)*M2^2))^-0.5;  % ??????
V = (1+2/((gamma-1)*M^2))^-0.5;
Vr = V*cos(beta);
Vt = -V*cos(beta)*tan(beta-delta);
y0 = [Vr;Vt];
t = [beta 0];
options = odeset('Events',@detectConeSurface,'AbsTol',1e-10);
[T,~] = ode45(@TEeq,t,y0,options,gamma);
tc = T(end);
end

function dy = TEeq(t,y,gamma)
% Taylor-Maccoll equation
% t = theta;
% y1 = Vr;
% y2 = Vtheta = dy1/dt;
% dy2/dt = num/den;
% num = -(gamma-1)/2*(1-y1^2-y2^2)*(2*y1+y2*cot(t))+y1*y2^2
% den = (gamma-1)/2*(1-y1^2-y2^2)-y2^2

tmp = (gamma-1)/2*(1-y(1)^2-y(2)^2);
num = -tmp*(2*y(1)+y(2)*cot(t))+y(1)*y(2)^2;
den = tmp-y(2)^2;
dy = [y(2);
      num/den];
end

function [value,isterminal,direction] = detectConeSurface(~,y,~)
value = y(2);
isterminal = 1;
direction = 0;
end

function delta = obliqueShock(gamma,M,beta)
delta = atan(2*cot(beta)*((M*sin(beta))^2-1)/(M^2*(gamma+cos(2*beta))+2));
if beta<=asin(1/M)
    delta = 0;
end
end