clear;
clc;

%% Power Flow calculation
Y=Ybus_new(case9_new_Sauer); % 9 bus system data obtained from MATPOWER
%Y=Ybus_new(case5_Overbye); 

%Y=Ybus_new(case14); 
result=runpf(case9_new_Sauer); % run ac power flow, in this case default NR is used
%result= runpf(case5_Overbye); 
%result=runpf(case14);
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
Xd=[0.06080; 0.11980; 0.18130];
R=[0;0;0];
%H=[23.64; 6.4; 3.01];
H=[13.64; 6.4; 3.01]; 
D=[0.0255; 0.00663; 0.00265]*10; 
%D=[9.6; 2.5; 1]; 
f0=60; 
w_syn=2*pi*f0; 
M=2*H/w_syn; 
gen_bus=result.gen(:, 1); 

%% machine data for 14 bus system 
% % Machine data 
% H=[5.1498; 6.54; 6.54; 5.06; 5.06];
% Xd=[0.2995; 0.185; 0.185; 0.232; 0.232];
% R=zeros(length(Xd), 1); 
% D=[2; 2; 2; 2; 2]*1; 
% f0=60; 
% w_syn=2*pi*f0; 
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
F=[4 6];
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
t_step=0.0005; 
t_SW=1; 
t_FC=1.0666; 
t_end=10; 
t=(0:t_end/t_step)*t_step; 

Ig=conj(Sg./V(1:length(result.gen(:, 1))));
E0=V(gen_bus)+Ig.*(R+1j*Xd);  % Machine terminal voltage
E_abs=abs(E0);

I0=Ybf*E0;
delta0=angle(E0)*180/pi;
w0=zeros(length(Xd), 1);
X_0=[angle(E0) w0] ; 

% Initialize power injection 
PG0=real(E0.*conj(I0)); 
PM=PG0;
QG0=imag(E0.*conj(I0));

YBUS(:, :, 1)=Ybf; 
YBUS(:, :, 2)=Ydf; 
YBUS(:, :, 3)=Yaf;
n_Gen=length(gen_bus); 
s=length(result.bus(:, 1)); 
n=n_Gen; 

%% Actual states using RK4 method 
kSW= t_SW/t_step+1;  
kend= t_end/t_step+1; 

% state 
x=zeros(2*n, kend); 
w=zeros(n, kend); 
w(:, 1)=w0; 
delta= zeros(n, kend); 
delta(:, 1)=delta0; 

% Measurements 
PG1=zeros(n, kend); 
QG1=zeros(n, kend); 
Vmag=zeros(s, kend); 
Vangle=zeros(s, kend); 
z=zeros(2*n+2*s, kend); 

for k=2:kend
    
    % Ybus and reconstruction matrix accodring to the requirement
    if k<=kSW
        ps=1;
    else 
        ps=3;
    end  
    Ybusm = YBUS(:,:,ps);
    RVm=RV(:, :, ps);
    
    % Generator voltage and current 
    E1=E_abs.*exp(1j*delta(:, k-1)); 
    I1=Ybusm*E1; 
    
    % Measurements: PG, QG
    PG1(:, k-1)=real(E1.*conj(I1)); 
    QG1(:, k-1)=imag(E1.*conj(I1));
    Vmag(:, k-1)=abs(RVm*E1); 
    Vangle(:, k-1)=angle(RVm*E1);
    
    % states using RK4 method 
    % k1
    k1_w=t_step*(M.^-1.*(PM-PG1(:, k-1)-(D.*((w(:, k-1)))))); 
    k1_delta=t_step*(w(:, k-1)); 
    
    
    % k2
    E2=E_abs.*exp(1j*(delta(:, k-1)+k1_delta/2)); 
    I2=Ybusm*E2; 
    PG2=real(E2.*conj(I2)); 
    k2_w= t_step*(M.^-1.*(PM-PG2-(D.*((w(:, k-1)+k1_w/2))))); 
    k2_delta=t_step*(w(:, k-1)+k1_w/2); 
    
    % k3
    E3=E_abs.*exp(1j*(delta(:, k-1)+k2_delta/2)); 
    I3=Ybusm*E3; 
    PG3=real(E3.*conj(I3)); 
    k3_w=t_step*(M.^-1.*(PM-PG3-(D.*((w(:, k-1)+k2_w/2))))); 
    k3_delta= t_step*(w(:, k-1)+k2_w/2); 
    
    % k4 
    E4=E_abs.*exp(1j*(delta(:, k-1)+k3_delta)); 
    I4=Ybusm*E4; 
    PG4=real(E4.*conj(I4)); 
    k4_w=t_step*(M.^-1.*(PM-PG4-(D.*((w(:, k-1)+k3_w))))); 
    k4_delta=t_step*(w(:, k-1)+k3_w); 
    
    w(:, k)=w(:, k-1)+(1/6)*(k1_w+2*k2_w+2*k3_w+k4_w); 
    delta(:, k)=delta(:, k-1)+(1/6)*(k1_delta+2*k2_delta+2*k3_delta+k4_delta);    
end 

% Actual state x, w in rad/sec, delta in rad 
X=[delta;  w_syn+w]; 
X=X(:, 1:kend-1); % because the last value was not able to update 

t=t(1:kend-1); 

for i=1:n
figure(i)%w1 in pu
plot(t, X(i, :), 'linewidth',2); 
% hold on 
% plot(t, X_est(:, i+n)+2*pi*f0, 'linestyle', '--', 'color', 'r', 'linewidth', 2); 
grid on 
xlabel('time(s)', 'fontsize', 15); 
ylabel('\delta in rad', 'fontsize', 15); 
%legend('\omega ', '\omega^e^s^t'); 
end 

% \delta with respect to first generator
figure(4)
subplot(2,1,1)
plot(t,X(2, :)-X(1, :), 'linewidth', 1.5)
grid on
ylabel('\delta_{2-1} [radian]', 'fontsize', 12)
title('Absolute coordinates', 'fontsize', 12)


subplot(2,1,2)
plot(t,X(3, :)-X(1, :), 'linewidth', 1.5)
grid on
ylabel('\delta_{3-1} [radian]', 'fontsize', 12)
xlabel('Time [s]', 'fontsize', 12)


% time=[t_SW; t_SW]; % fault occurance and clearance time, while determine DSE we don't care fault clearing time
% 
% options = odeset('RelTol',1e-9,'AbsTol',ones(2*n_Gen,1)*1e-9);
% [T X] = ode45(@(t,x) dynamic_system(t,x,M,D,YBUS,E_abs,PM,time,n_Gen),[t],X_0, options);
% 
% % %% Graph results
% % figure(1)
% % subplot(2,1,1)
% % plot(T,X(:,1:n_Gen)*180/pi, 'linewidth', 1.5)
% % grid on
% % ylabel('\delta_{G} [degree]', 'fontsize', 12)
% % title('Absolute coordinates', 'fontsize', 12)
% % subplot(2,1,2)
% % plot(T,X(:,n_Gen+1:2*n_Gen)+2*pi*f0, 'linewidth', 1.5)
% % grid on
% % ylabel('\omega_{s} [rad/sec]', 'fontsize', 12)
% % xlabel('Time [s]', 'fontsize', 12)
% 
% 
% % 
% % % results with respect to the center of inertia
% % MT = sum(M);
% % for k = 1:length(X(:,1))
% %     d_o = sum(X(k,1:n_Gen).*M')/MT;
% %     w_o = sum(X(k,n_Gen+1:2*n_Gen).*M')/MT;
% %     Xcoi(k,:) = X(k,:) - [d_o*ones(1,n_Gen),w_o*ones(1,n_Gen)];
% % end
% % 
% % figure(2)
% % subplot(2,1,1)
% % plot(T,Xcoi(:,1:n_Gen)*180/pi, 'linewidth', 1.5)
% % grid on
% % ylabel('\delta_{G} [degree]', 'fontsize', 15)
% % title('Center of Inertia Coordinate', 'fontsize', 15)
% % subplot(2,1,2)
% % plot(T,Xcoi(:,n_Gen+1:2*n_Gen)+2*pi*f0, 'linewidth', 1.5)
% % grid on
% % ylabel('\omega_{s} [rad/sec]', 'fontsize', 15)
% % xlabel('Time [s]', 'fontsize', 15)
% 
% % %% Graph results
% % figure(3)
% % subplot(2,1,1)
% % plot(T,X(:,2)-X(:,1), 'linewidth', 1.5)
% % grid on
% % ylabel('\delta_{2-1} [radian]', 'fontsize', 12)
% % title('Absolute coordinates', 'fontsize', 12)
% % subplot(2,1,2)
% % plot(T,X(:,3)-X(:,1), 'linewidth', 1.5)
% % grid on
% % ylabel('\delta_{3-1} [radian]', 'fontsize', 12)
% % xlabel('Time [s]', 'fontsize', 12)
% 
% 
% ps=1; 
% for k=1:length(T)
%     
%     if T(k)<=t_SW
%         ps=1;
% %     elseif (t_SW<T(k))&&(T(k)<t_FC)
% %         ps=2;
%     else 
%         ps=3; 
%     end  
%     Ybusm = YBUS(:,:,ps);
%     E1 = E_abs.*exp(1j*X(k, (1:n_Gen))');
%     Ibus = Ybusm*E1;
%     % measurements : PG, QG
%     PG1(:, k)=real(E1.*conj(Ibus)); 
%     QG1(:, k)=imag(E1.*conj(Ibus));
%     Vmag(:, k)=abs(RV(:, :, ps)*E1); 
%     Vangle(:, k)=angle(RV(:, :, ps)*E1);
%     z(:, k)=[PG1(:, k); QG1(:, k); Vmag(:, k); Vangle(:, k)]; 
%     
% end 
% 
% 
% 
% 
% % Actual state:
% kSW=t_SW/t_step+1; 
% kend=t_end/t_step+1; 
% 
% 
% % %Actual state x
% % %w in pu
% % %delta in rad
% 
% % t=t(1:kend-1); 
% % x=[1+w/w_syn; delta]; 
% % x=x(:, 1:kend-1); 
% 
% 
% 
% 
% %% Estimated State: 
% % Number of states and measurements 
% n=n_Gen;
% ns=2*n; 
% nm=2*n+2*s; 
% 
% % Covariance Matrix
% sig=1e-2; 
% P=sig^2*eye(ns);  % Error covariance matrix 
% Q=sig^2*eye(ns); % system noise covariance matrix 
% R=sig^2*eye(nm); % measurment noise covariance matrix 
% 
% kend=length(T); 
% % estimated state
% X_hat=zeros(2*n, kend); 
% 
% 
% % sigma points 
% X_sigma=zeros(ns, 2*ns); 
% 
% % constant values 
% PM_rep=repmat(PM, 1, 2*ns); 
% M_rep=repmat(M, 1, 2*ns); 
% D_rep=repmat(D, 1, 2*ns); 
% 
% %Unscented Kalman Filter (UKF) ALgorithm 
% for k=2:kend
%     
%     % Ybus and reconstruction matrix accodring to the requirement
%     if T(k-1)<=t_SW
%         ps=1;
%     else 
%         ps=3;
%     end  
%     Ybusm = YBUS(:,:,ps);
%     RVm=RV(:, :, ps);
%     
% 
%     % Sigma points
%     A=chol(ns*P); 
%     x_tilda=[A, -A]; 
%     X_sigma=repmat(X_hat(:, k-1), 1, 2*ns)+x_tilda; 
%     w_sigma=X_sigma(1:n, :);
%     delta_sigma=X_sigma(n+1:ns, :);
%     
% 
%     % Update sigma points 
%     E1=repmat(E_abs, 1, 2*ns).*exp(1j*delta_sigma); 
%     I1=Ybusm*E1; 
%     Pi1=real(E1.*conj(I1)); 
%     
% 
%     % k1
%     k1_w=t_step*(M_rep.^-1.*(PM_rep-Pi1-(D_rep.*(w_sigma./w_syn)))); 
%     k1_delta=t_step*w_sigma; 
% 
%     % k2
%     E2=repmat(E_abs, 1, 2*ns).*exp(1j*(delta_sigma+k1_delta/2)); 
%     I2=Ybusm*E2; 
%     Pi2=real(E2.*conj(I2)); 
%     k2_w=t_step*(M_rep.^-1.*(PM_rep-Pi2-(D_rep.*((w_sigma+k1_w/2)./w_syn)))); 
%     k2_delta=t_step*(w_sigma+k1_w/2); 
% 
%     % k3
%     E3=repmat(E_abs, 1, 2*ns).*exp(1j*(delta_sigma+k2_delta/2)); 
%     I3=Ybusm*E3; 
%     Pi3=real(E3.*conj(I3)); 
%     k3_w=t_step*(M_rep.^-1.*(PM_rep-Pi3-(D_rep.*((w_sigma+k2_w/2)./w_syn)))); 
%     k3_delta=t_step*(w_sigma+k2_w/2); 
% 
%     % k4
%     E4=repmat(E_abs, 1, 2*ns).*exp(1j*(delta_sigma+k3_delta)); 
%     I4=Ybusm*E4; 
%     Pi4=real(E4.*conj(I4)); 
%     k4_w=t_step*(M_rep.^-1.*(PM_rep-Pi4-(D_rep.*((w_sigma+k3_w)./w_syn)))); 
%     k4_delta=t_step*(w_sigma+k3_w); 
% 
%     w_sigma_up=w_sigma+(1/6)*(k1_w+2*k2_w+2*k3_w+k4_w); 
%     delta_sigma_up=delta_sigma+(1/6)*(k1_delta+2*k2_delta+2*k3_delta+k4_delta); 
%     w_sigma=w_sigma_up; 
%     delta_sigma=delta_sigma_up; 
%     X_sigma=[w_sigma; delta_sigma]; 
% 
%     % Priori State Estimate 
%     x_minus=(1/(2*ns))*(sum(X_sigma'))';
% 
%     % priori Covariance Matrix 
%     x_minus_rep=repmat(x_minus, 1, 2*ns); 
%     P_minus=(1/(2*ns))*(X_sigma-x_minus_rep)*(X_sigma-x_minus_rep)'+Q; 
% 
%     % New sigma points
%     A=chol(ns*P_minus); 
%     x_tilda=[A', -A']; 
%     X_sigma=repmat(X_hat(:, k-1), 1, 2*ns)+x_tilda; 
%     w_sigma=X_sigma(1:n, :); 
%     delta_sigma=X_sigma(n+1:ns, :); 
% 
%     % y sigma points
%     E=repmat(E_abs, 1, 2*ns).*exp(1j*delta_sigma); 
%     I=Ybusm*E; 
%     y1=real(E.*conj(I)); 
%     y2=imag(E.*conj(I)); 
%     y3=abs(RVm*E); 
%     y4=angle(RVm*E); 
%     y_sigma=[y1; y2; y3; y4]; 
% 
%     % y predict 
%     y_predict=(1/(2*ns))*(sum(y_sigma'))'; 
% 
%     % Covariance of predicted measurements Py
%     y_predict_rep=repmat(y_predict, 1, 2*ns); 
%     P_y=(2*ns)^-1*(y_sigma-y_predict_rep)*(y_sigma-y_predict_rep)'+R; 
% 
%     % Cross Covariance Pxy 
%     P_xy=(2*ns)^-1*(X_sigma-x_minus_rep)*(y_sigma-y_predict_rep)'; 
% 
%     % Measurement update of state estimate 
%     K=P_xy*P_y^-1; 
%     y_PG=PG1(:, k); 
%     y_QG=QG1(:, k); 
%     y_Vmag=Vmag(:, k); 
%     y_Vangle=Vangle(:, k); 
%     v=sig^2*rand(nm, 1); 
%     y=[y_PG; y_QG; y_Vmag; y_Vangle]+v; 
%     X_hat(:, k)=x_minus+K*(y-y_predict); 
%     P=P_minus-K*P_y*K';        
% 
%    
% end 
% 
% % Estimated State 
% % w in pu 
% % delta in rad 
% X_est=[ X_hat(n+1:ns, :); X_hat(1:n, :);]; 
% X_est=transpose(X_est(:, 1:kend)); 
% 
% %% Error calculation  
% %Overall Estimation Error 
% err_est=(1/size(X, 2))*sum((1/ns)*sum(abs(X-X_est))); 
% 
% % Estimation Error for w
% err_est_w=(1/size(X, 2))*sum((1/n)*sum(abs(X(1:n, :)-X_est(1:n, :)))); 
% 
% % Estimated error for delta 
% err_est_delta=(1/size(X, 2))*sum((1/n)*sum(abs(X(n+1:2*n, :)-X_est(n+1:2*n, :)))); 
% 
% %% Plot 
% t=t(1:kend); 
% % Fifuewa:
% for i=1:n
% figure(i)%w1 in pu
% plot(t, X(:, i), 'linewidth',2); 
% hold on 
% plot(t, X_est(:, i), 'linestyle', '--', 'color', 'r', 'linewidth', 2); 
% grid on 
% xlabel('time(s)', 'fontsize', 15); 
% ylabel('\delta in pu', 'fontsize', 15); 
% legend('\delta', '\delta^e^s^t'); 
% end 
% 
% 
% 
% for i=1:n
% figure(n+i)%w1 in pu
% plot(t, X(:, i+n)+2*pi*f0, 'linewidth',2); 
% hold on 
% plot(t, X_est(:, i+n)+2*pi*f0, 'linestyle', '--', 'color', 'r', 'linewidth', 2); 
% grid on 
% xlabel('time(s)', 'fontsize', 15); 
% ylabel('\omega in pu', 'fontsize', 15); 
% legend('\omega ', '\omega^e^s^t'); 
% end 
% 
% 
% % \delta with respect to first generator
% figure(7)
% subplot(2,1,1)
% plot(t,X(:,2)-X(:,1), 'linewidth', 1.5)
% hold on 
% plot(t, X_est(:, 2)- X_est(:, 1), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
% grid on
% ylabel('\delta_{2-1} [radian]', 'fontsize', 12)
% title('Absolute coordinates', 'fontsize', 12)
% legend('\delta_{2-1} ', '\delta_{2-1}^e^s^t'); 
% 
% subplot(2,1,2)
% plot(t,X(:,3)-X(:,1), 'linewidth', 1.5)
% hold on 
% plot(t, X_est(:, 3)-X_est(:, 1), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
% grid on
% ylabel('\delta_{3-1} [radian]', 'fontsize', 12)
% xlabel('Time [s]', 'fontsize', 12)
% legend('\delta_{3-1} ', '\delta_{3-1}^e^s^t'); 




