%% SP for COMM HW 1 Q 1
% Instantaneous MIMO model

% basic params
c = 3e8;         
fc = 1e9;      % centre frequency, Hz
lambda = c/fc;   % wavelength, m

%% Question 1A

dist = 0.15;     % 15 cm
Delta = dist/lambda;

M1 = 5;
theta1 = 45;
a = gen_a(M1, Delta, theta1);
disp(a)

%% Question 1B
M2 = 3;
w = ones(M2,1);
theta_range = [-90 90];
[y, thetas] = spat_response(w,M2,Delta, theta_range);

%% Question 1C
M3 = 5;
theta3 = [0 30];    % deg
SNR = 3;    % dB
N = 20;

X = gen_data(M3,N,Delta,theta3,SNR);

%% Question 2

% make Spatial Response plots
% 1) One source, fixed w and Delta, var M
M_vec = [2 3 7];
figure(1);
hold on;
for k = 1:3
M0 = M_vec(k);
w0 = ones(M0,1);
[y, thetas] = spat_response(w0,M0,Delta, theta_range);
plot(thetas, abs(y));
end
hold off;
title("Spatial response for fixed w and Delta="+Delta);
xlabel('theta [deg]');
legend('M=2','M=3','M=7');

% 2) One source, fixed w and M, var Delta
Delta_vec = [0.5 1 2];
figure(2);
hold on;
for k = 1:3
Delta0 = Delta_vec(k);
M = 7;
w = ones(M,1);
[y, thetas] = spat_response(w,M,Delta0, theta_range);
plot(thetas, abs(y));
end
hold off;
title("Spatial response for fixed w and M="+M);
xlabel('theta [deg]');
legend('Delta=0.5','Delta=1','Delta=2');


% 3) One source, fixed w and Delta, var M
%    not at 0 deg
M_vec = [2 3 7];
theta_source = 30;
figure(3);
hold on;
for k = 1:3
M0 = M_vec(k);
w0 = gen_a(M0,Delta,theta_source);
[y, thetas] = spat_response(w0,M0,Delta, theta_range);
plot(thetas, abs(y));
end
hold off;
title("Spatial response for source at 30 deg and Delta="+Delta);
xlabel('theta [deg]');
legend('M=2','M=3','M=7');


% 3) Two sources, fixed w and Delta and M, var theta
M = 3;
Delta = 0.5;
theta_arr = {[0 10],[0 30]};
for k = 1:2
figure(3+k);
hold on;
theta_sources = theta_arr{k};
theta1 = theta_sources(1);
theta2 = theta_sources(2);

w1 = gen_a(M, Delta, theta1);
w2 = gen_a(M, Delta, theta2);

w0 = [w1 w2]; 

[y, thetas] = spat_response(w0,M,Delta,theta_range);
plot(thetas, abs(y));
hold off;
title("Spatial response for fixed M="+M+", Delta="+Delta+"for two sources");
xlabel('theta [deg]');
str = {'Angle of sources (theta)',num2str(theta_sources)};
annotation('textbox',[0.15 0.5 0.3 0.3],'String',str,'FitBoxToText','on');
end

%% Question 3


Delta5 = 0.5;    % antenna separation
M5 = 5;          % antennas
N5 = 10;         % samples
theta5 = [30]; % deg
SNR5 = 3;        % dB
count = 0;
X_svd = {};

% varying M
M_var = [3, 7, 12];
for i = 1:3
    M5_var = M_var(i);
    X = gen_data(M5_var,N5,Delta5,theta5,SNR5);
    X_svd{count+i} = svd(X);  
end
count = 3;

% varying theta for 2 sources
theta_var = {[0 10],[0 30],[0 45]};
for i = 1:3
    theta5_var = theta_var{i};
    X = gen_data(M5,N5,Delta5,theta5_var,SNR5);
    X_svd{count+i} = svd(X);
end
count = 6;

% varying number of sources
theta_var = {[0 10],[0 10 20],[0 10 20 30]};
for i = 1:3
    theta5_var = theta_var{i};
    X = gen_data(M5,N5,Delta5,theta5_var,SNR5);
    X_svd{count+i} = svd(X);
end
count = 9;

% varying number of samples
N_var = [2,3,4];
for i = 1:3
    N5_var = N_var(i);
    X = gen_data(M5,N5,Delta5,theta5_var,SNR5);
    X_svd{count+i} = svd(X);
end
count = 12;

% varying SNR
SNR_var = [3 10 20];
for i = 1:3
    SNR5_var = SNR_var(i);
    X = gen_data(M5,N5,Delta5,theta5,SNR5_var);
    X_svd{count+i} = svd(X);
end

figure(6);
str_tit = {'(a) varying M','(b) varying theta, 2 sources','(c) varying number of sources','(d) varying N','(e) varying SNR'};
count = 0;
for h = 1:5
   AX(i) = nexttile;
   hold on;
   p1 = plot(X_svd{count+1},'k-o');
   p2 = plot(X_svd{count+2},'b-o');
   p3 = plot(X_svd{count+3},'r-o');
   count = count+3;
   title(str_tit{h});
   legend('show');
   hold off;
end


%% Functions

function a = gen_a(M,Delta,theta)
    % Generate uniform linear array response vector for model
    %
    % INPUT:
    % M : number of elements
    % Delta : spacing between elements
    % theta : direction of incident wave in DEGREES
    %
    % OUTPUT:
    % a : array response
    

    ang = deg2rad(theta);   % convert to degree
    a = zeros(M,1);         % preallocate

    phi = exp(2j*pi*Delta*sin(ang)); 

    for i=1:M
        a(i) = phi^(i-1);   % array response vector
    end
    
    % FUNCTION CHECK:
    % elementPos = (0:dist:dist*(M-1));
    % sv = steervec(elementPos/lambda,theta);
    % disp(sv)
end

function [y,thetas] = spat_response(w, M, Delta, theta_range)
% generate spatial response y =  |w*A(theta)|
% INPUT:
% M : number of antennas
% w : given beamformer
% Delta : antenna spacing (electrical)
% theta_range : directions range [min max]
%
% OUTPUT:
% y : spatial response

theta_res = 1;
thetas = (theta_range(1):theta_res:theta_range(2));

for i=1:length(thetas)
    A(:,i) = gen_a(M, Delta, thetas(i));
end

y = w'*A;
end

function X = gen_data(M,N,Delta,thetas,SNR)
    % generate data matrix X = AS + N
    % INPUT:
    % M : number of antennas
    % N : number of samples
    % Delta : antenna spacing (electrical)
    % thetas : directions (vector)
    % SNR : in dB
    %
    % OUTPUT:
    % X : data matrix

    % generate array response matrix
    d = length(thetas);
    for i=1:d
        A(:,i) = gen_a(M, Delta, thetas(i));
    end

    % convert SNR to ration
    SNR_W = 10^(SNR/10);

    % generate zero-mean random Gaussian complex source data
    Source = randn([d,N]) + 1j*randn([d,N]);
    Source = Source - mean(Source,'all');

    % generate zero-mean random complex Gaussian noise
    Noise = randn([M,N]) + 1j*randn([M,N]);
    Noise = Noise - mean(Noise,'all');

    var_noise = var(A*Source,0,'all')/SNR_W;
    std_noise = sqrt(var_noise);
    scale_noise = std(Noise,0,'all')/std_noise; 
    Noise = Noise/scale_noise;

    X = A*Source + Noise;
end


