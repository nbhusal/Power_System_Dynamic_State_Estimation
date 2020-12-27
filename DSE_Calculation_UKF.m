clear;
clc;

%% Power Flow calculation
% Y=Ybus_new(case9_new_Sauer); % 9 bus system data obtained from MATPOWER
% result=runpf(case9_new_Sauer); % run ac power flow, in this case default NR is used

% Y=Ybus_new(case5_Overbye); 
% result= runpf(case5_Overbye); 

% Y=Ybus_new(case14); 
% result=runpf(case14);

Y=Ybus_new(case39);
result=runpf(case39);

% Y=Ybus_new(case145);
% result=runpf(case145);

Vmag=result.bus(:, 8); % Pu voltage magnitude of each buses 
Vph=result.bus(:, 9); % angle in degree
V=Vmag.*exp(1j*Vph*pi/180); 
P_jQ=conj(V).*(Y*V); % Net Power at each node
S=conj(P_jQ);
S=S/100; 
Sg=result.gen(:, 2)+1j*result.gen(:, 3); 
Sg=Sg/100;


%% machine data for 9 bus system
% Xd=[0.06080; 0.11980; 0.18130];
% R=[0;0;0];
% H=[23.64; 6.4; 3.010];
% M=H/(pi*60); 
%D=[0.0125;0.0034;0.0016];
% 
%% Data of 9 bus system from Peter Sauer.
% Xd=[0.06080; 0.11980; 0.18130];
% R=[0;0;0];
% H=[23.64; 6.4; 3.01];
% %H=[13.64; 6.4; 3.01]; 
% D=[0.0255; 0.00663; 0.00265]; 
% %D=[9.6; 2.5; 1]; % If we use this value need to devide the D term by 2*pi*60 
% f0=60; 
% w_syn=2*pi*f0; 
% M=2*H/w_syn; 
% gen_bus=result.gen(:, 1); 

%% machine data for 14 bus system 
% % Machine data 
% H=[5.1498; 6.54; 6.54; 5.06; 5.06];
% Xd=[0.2995; 0.185; 0.185; 0.232; 0.232];
% R=zeros(length(Xd), 1); 
%  
% f0=60; 
% w_syn=2*pi*f0; 
% D=[2; 2; 2; 2; 2]/w_syn;
% 
% M=2*H/w_syn; 
% 
% gen_bus=result.gen(:, 1); 

%% Overbye data for 5 bus system 
% Xd=[0.05; 0.025]; 
% R=[0; 0]; 
% H=[]; 
% D=[]; 
% f0=60; 
% w_syn=2*pi*f0; 
% M=2*H/w_syn; 
% gen_bus=result.gen(:, 1); 

%% Case 39 bus data 
Xd=[0.006; 0.0697; 0.0531; 0.0436; 0.132; 0.05; 0.049; 0.057; 0.057; 0.031]; 
H=[500; 30.3; 35.8;28.6; 26; 34.8; 26.4; 24.3; 34.5; 42]; 
R=zeros(length(Xd), 1); 
f0=60; 
w_syn=2*pi*f0; 

D=[0; 0;0 ;0; 0; 0; 0; 0; 0; 0]; 
D=D/w_syn;

M=2*H/w_syn; 

gen_bus=result.gen(:, 1); 

%% case 145



%% calculate Y22
Y22=diag(1./(1j*Xd)); 

%% Calculation of Y11
SL=result.bus(:, 3)+1j*result.bus(:, 4); 
SL=SL/100; 
YL=conj(SL)./(abs(V).^2); % 
Y11=Y+diag(YL);
Y11(gen_bus, gen_bus)=Y11(gen_bus, gen_bus)+Y22;
% 

%% Calculation of Y12 and Y21
% Calculation of Y12


Y12=zeros(length(result.bus(:,1)), length(result.gen(:,1)));
%Y12(gen_bus, gen_bus)=Y12(gen_bus, gen_bus)-Y22; 

for i=1:length(result.bus(:,1))
    for k=1:length(result.gen(:,1))
        q=result.gen(k,1);
        if i==q
            Y12(q,k)=-1/(R(k)+Xd(k)*1j);
        end 
    end 
end 

Y21=transpose(Y12);
%% Calculation of reduced matrix before fault
Ybf=Y22-Y21*inv(Y11)*Y12 ;

% Bus Reconstruction matrix 
RV(:, :, 1)=-inv(Y11)*Y12;


%% Enter fault here to calculate the afterfault and during fault reduced
% matrices
f11=4;
F=[4 14];
f1=F(1);
f2=F(2);
%% during fault
Y11df=Y11; 
Y11df(f11, :)=[];
Y11df(:,f11)=[];
Y12df=Y12; 
Y12df(f11, :)=[];
Y21df=transpose(Y12df);
% during fault reduced matrics
Ydf=Y22-Y21df*inv(Y11df)*Y12df;

RV(:, :, 2)=zeros(size(RV(:, :, 1)));
RV(1:end-1, :, 2)=RV(1:end-1, :, 2)-inv(Y11df)*Y12df;


%% afterfault Y11
 Y11after=Y11;
 Y11after(f1,f2)=0;
 Y11after(f2,f1)=0;
for i=1:length(result.branch(:,1)) 
    if (f1==result.branch(i,1)&& f2==result.branch(i,2))||(f2==result.branch(i,1)&& f1==result.branch(i,2))
        Y11after(f1,f1)=Y11after(f1,f1)-result.branch(i,5)*1j/2-1/(result.branch(i,3)+result.branch(i,4)*1j);
        Y11after(f2,f2)=Y11after(f2,f2)-result.branch(i,5)*1j/2-1/(result.branch(i,3)+result.branch(i,4)*1j);
    end
end 

% Afterfault reduced matrix is 
Yaf=Y22-Y21*inv(Y11after)*Y12 ;
%RV_af=-inv(Y11after)*Y12 ; 

RV(:, :, 3)=-inv(Y11after)*Y12; 
%% Initialization
deltt=0.0005; 
t_SW=1; 
t_FC=1.0333; 
t_max=10; 


Ig=conj(Sg./V(1:length(result.gen(:, 1))));
E0=V(gen_bus)+Ig.*(R+1j*Xd);  % Machine terminal voltage
E_abs=abs(E0);

I0=Ybf*E0;
delta0=angle(E0)*180/pi;
w0=zeros(length(Xd), 1);
X_0=[angle(E0); w0]; 

% Initialize power injection 
PG0=real(E0.*conj(I0)); 
PM=PG0;
QG0=imag(E0.*conj(I0));

YBUS(:, :, 1)=Ybf; 
YBUS(:, :, 2)=Ydf; 
YBUS(:, :, 3)=Yaf;
n=length(gen_bus); 
s=length(result.bus(:, 1)); 




%% Estimated State: 
% Number of states and measurements 
ns=2*n; 
nm=2*n+2*s; 



% Covariance Matrix
sig=1e-2; 
P=sig^2*eye(ns);  % Error covariance matrix 
Q=sig^2*eye(ns); % system noise covariance matrix 
R=sig^2*eye(nm); % measurment noise covariance matrix 

X_hat=X_0;
X_est=[]; 
X_mes=[]; % Initial statel 

% constant values  
W=ones(ns*2,1)/(2*ns);
RMSE=[]; 
%Unscented Kalman Filter (UKF) ALgorithm 
for k=0:deltt:t_max
    % Ybus and reconstruction matrix accodring to the requirement
    if k<t_SW
        ps=1;
    elseif (t_SW<k)&&(k<=t_FC)
        ps=2;  
    else 
        ps=3; 
    end  
    Ybusm = YBUS(:,:,ps);
    RVm=RV(:, :, ps);
    
    [~, X] = ode45(@(t,x) dynamic_system(t,x,M,D,Ybusm,E_abs,PM,n),[k k+deltt],X_0);
    X_0=transpose(X(end, :));
    X_mes=[X_mes X_0];
    
    %determine the measurements 
    E1=E_abs.*exp(1j*X_0(1:n)); 
    I1=Ybusm*E1; 
    PG=real(E1.*conj(I1)); 
    QG=imag(E1.*conj(I1)); 
    Vmag=abs(RVm*E1); 
    Vangle=angle(RVm*E1); 
    z=[PG; QG; Vmag; Vangle]; 
      
    
    % Sigma points for X
    root=chol(ns*P); 
    X_tilde=[root, -root]; 
    X_sigma=repmat(X_hat, 1, 2*ns)+X_tilde;
    
    
     
    % update estimate    
    xbreve=RK4(n, deltt, E_abs, ns, X_sigma, PM, M, D, Ybusm); % Using RK4

%     xbreve=X_sigma;
%     for i=1:2*ns % using Ode
%         [T1, X1] = ode45(@(t,x) dynamic_system(t,x,M,D,Ybusm,E_abs,PM,n),[k k+deltt], X_sigma(:, i));
%         xbreve(:, i)=transpose(X1(end, :));
%     end
% %     
    X_hat=xbreve*W;

    % priori Covariance Matrix   
    x_hat_rep=repmat(X_hat, 1, 2*ns);
    P=(1/(2*ns))*(xbreve-x_hat_rep)*transpose((xbreve-x_hat_rep))+Q;
    
    % New sigma points
    root1=chol(ns*P, 'lower'); 
    
    X_tilde1=[root1, -root1];  
    X_sigma=repmat(X_hat, 1, 2*ns)+X_tilde1;
    
    % sigma points 
    E11=repmat(E_abs, 1, ns*2).*exp(1j*X_sigma(1:n, :)); 
    I11=Ybusm*E11;    
    PG11=real(E11.*conj(I11)); 
    QG11=imag(E11.*conj(I11)); 
    Vmag11=abs(RVm*E11); 
    Vangle11=angle(RVm*E11); 
    zbreve=[PG11; QG11; Vmag11; Vangle11]; 
    
    % predict zhat 
    zhat=zbreve*W;
    
    % Covariance of predicted measurements Pz and cross covariance Pxz
    zhat_rep=repmat(zhat, 1, 2*ns); 
    Pz=(1/(2*ns))*(zbreve-zhat_rep)*transpose((zbreve-zhat_rep))+R;
    
    Pxz=(1/(2*ns))*(X_sigma-x_hat_rep)*transpose((zbreve-zhat_rep));

   
    
    % Measurement update of state estimate 
    K=Pxz/Pz;
    X_hat=X_hat+K*(z-zhat); 
    P=P-K*Pz*transpose(K); 
    
    X_est=[X_est, X_hat]; 
     RMSE=[RMSE, sqrt(trace(P))];
end 

save('39_RMSE_UKF.mat', 'RMSE')


%% Plots
t= (0:deltt:t_max);
for i=1:n
figure(i)
subplot(2,1,1)
plot(t,X_mes(i, :), 'linewidth', 1.5)
hold on 
plot(t, X_est(i, :), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
grid on
ylabel(sprintf('Angle_{%d}', i), 'fontsize', 12)
xlabel('time(s)', 'fontsize', 15); 
title('Actual Vs Estimated \delta with UKF', 'fontsize', 12)
legend(sprintf('Angle_{%d, Actual} ',i), sprintf('Angle_{%d, UKF}', i), 'fontsize', 12, 'Location', 'northwest'); 

subplot(2,1,2)
plot(t,X_mes(i+n, :), 'linewidth', 1.5)
hold on 
plot(t, X_est(i+n, :), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
grid on
ylabel(sprintf('Speed_{%d}', i), 'fontsize', 12)
xlabel('time(s)', 'fontsize', 15); 
title('Actual Vs Estimated \omega with UKF', 'fontsize', 12)
legend(sprintf('Speed_{%d, Actual} ',i), sprintf('Speed_{%d, UKF}', i), 'fontsize', 12, 'Location', 'northwest'); 
end

% % results with respect to the center of inertia
% MT = sum(M);
% for k = 1:length(X_est(1, :))
%     d_oe = sum(X_est(1:n, k).*M')/MT;
%     d_oa = sum(X_mes(1:n, k).*M')/MT;
%     
%     w_oe = sum(X_est(n+1:2*n, k).*M')/MT;
%     w_oa = sum(X_mes(n+1:2*n, k).*M')/MT;
%     
%     Xcoie(k, :) = X_est(:, k) - [d_oe, w_oe]';
%     Xcoia(k, :) = X_mes(:, k) - [d_oa, w_oa]';
% 
% end

% %% Plots
% t= (0:deltt:t_max);
% for i=1:n
% figure(i+n)
% subplot(2,1,1)
% plot(t,Xcoia(:, i), 'linewidth', 1.5)
% hold on 
% plot(t, Xcoie(:, i), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
% grid on
% ylabel(sprintf('Angle_{%d}', i), 'fontsize', 12)
% xlabel('time(s)', 'fontsize', 15); 
% title('Measured Vs Eistimated \delta with respect to COI for UKF', 'fontsize', 12)
% legend(sprintf('Angle_{%d, Actual} ',i), sprintf('Angle_{%d, UKF}', i), 'fontsize', 12); 
% 
% subplot(2,1,2)
% plot(t,Xcoia(:, i+n), 'linewidth', 1.5)
% hold on 
% plot(t, Xcoie(:, i+n), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
% grid on
% ylabel(sprintf('Speed_{%d}', i), 'fontsize', 12)
% xlabel('time(s)', 'fontsize', 15); 
% title('Measured Vs Eistimated \omega with respect to COI for UKF', 'fontsize', 12)
% legend(sprintf('Speed_{%d, Actual} ',i), sprintf('Speed_{%d, UKF}', i), 'fontsize', 12); 
% end



% figure(i+n)
% subplot(2,1,1)
% plot(t,X_mes(3, :)-X_mes(1, :), 'linewidth', 1.5)
% hold on 
% plot(t, X_est(3, :)-X_est(1, :), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
% grid on
% ylabel(sprintf('Angle_{%d}', i), 'fontsize', 12)
% xlabel('time(s)', 'fontsize', 15); 
% title('Measured Vs Eistimated \delta with UKF', 'fontsize', 12)
% legend(sprintf('Angle_{%d, Actual} ',i), sprintf('Angle_{%d, UKF}', i), 'fontsize', 12); 
% 
% subplot(2,1,2)
% plot(t,X_mes(2, :)-X_mes(1, :), 'linewidth', 1.5)
% hold on 
% plot(t, X_est(2, :)-X_est(1, :), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
% grid on
% ylabel(sprintf('Angle_{%d}', i), 'fontsize', 12)
% xlabel('time(s)', 'fontsize', 15); 
% title('Measured Vs Eistimated \delta with UKF', 'fontsize', 12)
% legend(sprintf('Angle_{%d, Actual} ',i), sprintf('Angle_{%d, UKF}', i), 'fontsize', 12);




%  USE ode45 to solve this differential equation 
function dx = dynamic_system(~,x,M,D, Ybusm,Vo,Pm,NumG)  % dynamics of the sin control system
Vg = Vo.*exp(1j*x(1:NumG));
Ibus = Ybusm*Vg;
S = conj(Ibus).*Vg;
Pe = real(S);
dx = zeros(2*NumG,1);
dx(1:NumG) = x(NumG+1:2*NumG);
dx(NumG+1:2*NumG) = (Pm-Pe)./M-D.*x(NumG+1:2*NumG)./M;
end 

function xbreve=RK4(n, deltt, E_abs, ns, X_sigma, PM, M, D, Ybusm)
    % update sigma point 
    E1=repmat(E_abs, 1, 2*ns).*exp(1j*X_sigma(1:n, :)); 
    I1=Ybusm*E1; 
    PG1=real(E1.*conj(I1)); 
    PM_rep=repmat(PM, 1, 2*ns);
    M_rep=repmat(M, 1, 2*ns);
    D_rep=repmat(D, 1, 2*ns); 
    
    % k1
    k1_w=deltt*(M_rep.^-1.*(PM_rep-PG1-(D_rep.*(X_sigma(n+1:ns, :)))));
    k1_delta=deltt*X_sigma(n+1:ns, :);
    
    % k2
    E2=repmat(E_abs, 1, 2*ns).*exp(1j*(X_sigma(1:n, :)+k1_delta/2)); 
    I2=Ybusm*E2; 
    PG2=real(E2.*conj(I2));     
    k2_w=deltt*(M_rep.^-1.*(PM_rep-PG2-(D_rep.*((X_sigma(n+1:ns, :)+k1_w/2)))));
    k2_delta=deltt*(X_sigma(n+1:ns, :)+k1_w/2); 
    
    % k3
    E3=repmat(E_abs, 1, 2*ns).*exp(1j*(X_sigma(1:n, :)+k2_delta/2)); 
    I3=Ybusm*E3; 
    PG3=real(E3.*conj(I3));     
    k3_w=deltt*(M_rep.^-1.*(PM_rep-PG3-(D_rep.*((X_sigma(n+1:ns, :)+k2_w/2)))));
    k3_delta=deltt*(X_sigma(n+1:ns, :)+k2_w/2); 
    
    % k2
    E4=repmat(E_abs, 1, 2*ns).*exp(1j*(X_sigma(1:n, :)+k3_delta)); 
    I4=Ybusm*E4; 
    PG4=real(E4.*conj(I4));     
    k4_w=deltt*(M_rep.^-1.*(PM_rep-PG4-(D_rep.*((X_sigma(n+1:ns, :)+k3_w)))));
    k4_delta=deltt*(X_sigma(n+1:ns, :)+k3_w); 
    
    xbreve(1:n, :)=X_sigma(1:n, :)+ (k1_delta+2*k2_delta+2*k3_delta+k4_delta)/6;
    xbreve(n+1:ns, :)=X_sigma(n+1:ns, :) + (k1_w+2*k2_w+2*k3_w+k4_w)/6;
end

