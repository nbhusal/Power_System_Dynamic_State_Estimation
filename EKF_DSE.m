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
H=[23.64; 6.4; 3.01];
%H=[13.64; 6.4; 3.01]; 
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
deltt=0.0005; 
t_SW=1; 
t_FC=1.0666; 
t_max=10; 


Ig=conj(Sg./V(1:length(result.gen(:, 1))));
E0=V(gen_bus)+Ig.*(R+1j*Xd);  % Machine terminal voltage
E_abs=abs(E0);

I0=Ybf*E0;
delta0=angle(E0)*180/pi;
w0=zeros(length(Xd), 1);
X_0=[angle(E0); w0] ; 

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
X_est=X_0; 
X_mes=X_0; % Initial statel 

% constant values 
PM_rep=repmat(PM, 1, 2*ns); 
M_rep=repmat(M, 1, 2*ns); 
D_rep=repmat(D, 1, 2*ns); 
W=ones(ns*2,1)/(2*ns);

%Unscented Kalman Filter (UKF) ALgorithm 
for k=0:deltt:t_max
    k
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
    
    [T, X] = ode45(@(t,x) dynamic_system(t,x,M,D,Ybusm,E_abs,PM,n),[0 deltt],X_0);
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
    A=chol(ns*P); 
    x_tilde=[A, -A]; 
    X_sigma=repmat(X_hat, 1, 2*ns)+x_tilde;
   
    xbreve=X_sigma;
    for i=1:2*ns
        [T1, X1] = ode45(@(t,x) dynamic_system(t,x,M,D,Ybusm,E_abs,PM,n),[0 deltt], X_sigma(:, i));
        xbreve(:, i)=transpose(X1(end, :));
    end
    
    xhat=xbreve*W; 

    % priori Covariance Matrix 
    x_minus_rep=repmat(xhat, 1, 2*ns); 
    P=(1/(2*ns))*(X_sigma-x_minus_rep)*(X_sigma-x_minus_rep)'+Q; 

    % New sigma points
    A=chol(ns*P); 
    x_tilde=[A, -A]; 
    X_sigma=repmat(X_hat, 1, 2*ns)+x_tilde; 
    
    for i=1:2*ns
    E11=E_abs.*exp(1j*X_sigma(1:n, i)); 
    I11=Ybusm*E11; 
    PG1=real(E11.*conj(I1)); 
    QG1=imag(E11.*conj(I1)); 
    Vmag1=abs(RVm*E11); 
    Vangle1=angle(RVm*E1); 
    zbreve(:, i)=[PG1; QG1; Vmag1; Vangle1]; 
    end
    
    zhat=zbreve*W;
    zhat_rep=repmat(zhat, 1, 2*ns);
    xhat_rep=repmat(xhat, 1, 2*ns);
    
    P_y=(1/(2*ns))*(zbreve-zhat_rep)*transpose((zbreve-zhat_rep))+R;
    % Cross Covariance Pxy 
    P_xy=(1/2*ns)*(X_sigma-x_minus_rep)*transpose((zbreve-zhat_rep)); 
    
    % Measurement update of state estimate 
    K=P_xy/P_y;
    X_hat=X_hat+K*(z-zhat); 
    P=P-K*P_y*transpose(K); 
    
    X_est=[X_est, X_hat]; 

   
end 



%% Plots
t= (0:deltt:t_max);
for i=1:n
figure(i)
subplot(2,1,1)
plot(t,X_mes(i, :), 'linewidth', 1.5)
hold on 
plot(t, X_est(i, :), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
grid on
ylabel(sprintf('X_{%d}', i), 'fontsize', 12)
xlabel('time(s)', 'fontsize', 15); 
title('Measured Vs Eistimated \delta with UKF', 'fontsize', 12)
legend(sprintf('X_{%d, mes} ',i), sprintf('X_{%d, est}', i), 'fontsize', 12); 

subplot(2,1,2)
plot(t,X_mes(i+n, :), 'linewidth', 1.5)
hold on 
plot(t, X_est(i+n, :), 'linestyle', '--', 'color', 'r', 'linewidth', 2);
grid on
ylabel(sprintf('X_{%d}', i+n), 'fontsize', 12)
xlabel('time(s)', 'fontsize', 15); 
title('Measured Vs Eistimated \omega with UKF', 'fontsize', 12)
legend(sprintf('X_{%d, mes} ',i+n), sprintf('X_{%d, est}', i+n), 'fontsize', 12); 
end



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
