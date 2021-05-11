%% SP for COMM HW 1 Q 2
% Convolutive Model

%% Question 1 (a)
L = 3;
syms t
t = 0:.1:4;

tres = .1;                          % time delay resolution
tau_array = (0:tres:1-tres);        % delay values
idx=randperm(length(tau_array),1);  % get random index
tau = tau_array(idx);               % get random delay 'tau'

% Alternatively, choose a specific delay value by commenting the tau above
% and entering a value below:
% tau = 0.5;    % choose from 0, 0.1, ... 0.8, 0.9

% choose arbitrary P, L values for now
P = 15;
L = 3;

% sample pulse with delay
g = pulse(tau,L,P);

% plot reference pulse and delayed pulse
t22 = linspace(0,3,L*P);        % new time axis
for k = 1
    tri = @(t) triangularPulse(0,2,t);
    figure(1);
    plot(t,tri(t));
    hold on;
    plot(t22,g);
    ylabel('g(t)'); xlabel('t');
    str = {'Time Delay tau',num2str(tau)};
    annotation('textbox',[0.2 0.5 0.3 0.3],'String',str,'FitBoxToText','on');
    title('reference pulse and delayed pulse');
end

%% Question 1(b)
r = 2;                                  % arbitrary number of paths
idxs = randperm(length(tau_array),r);   % get random index
tau_vec = tau_array(idxs);              % generate time delay vector

beta = 1;                               % gains, column vector

h = channel(tau_vec, beta, L, P);

%% Question 1 (c)
qpsk_symb = [1 1j -1 -1j];
N = 3;
s = source(N, qpsk_symb)';

%% Question 1 (d)
[x, xx] = gen_data(h,s,P,N);

%% Question 2
% given param values
r2 = 2;
tau2 = [0.1 0.6]';
phi = [0 15];       % arbitrary, deg
beta2 = [1*exp(1j*deg2rad(phi(1))) 0.7*exp(1j*deg2rad(phi(2)))]';
P2 = 5;
N2 = 50;
% derived values
t2 = linspace(0,N2-1,N2);
L2 = N2+1;
h2 = channel(tau2, beta2, L2, P2)'; % channel
s2 = source(N2, qpsk_symb)';        % source
[x2, xx2] = gen_data(h2,s2,P2,N2);  % data

% plot results
lenx = length(x2);
t_ax2 = linspace(0,lenx+1,lenx);
figure(2); plot(t_ax2, real(x2)); hold on; plot(t_ax2, imag(x2))
legend('real','imag');title('real and imaginary parts of x2')
xlim([0 70]);

lenh = length(h2);
t_ax3 = linspace(0,lenh+1,lenh);
figure(3); plot(t_ax3, real(h2)); hold on; plot(t_ax3, imag(h2))
legend('real','imag');title('real and imaginary parts of h2')
xlim([0 15]);

fprintf('-------------------\n');
fprintf('Source s2 length: %d\n', length(s2));
fprintf('Channel h2 length: %d\n', length(h2));
fprintf('Convolution output xx2 length: %d\n', length(xx2));
fprintf('Truncated data vector x2 length: %d\n', length(x2));

%% Question 3
%  reshape
X = reshape(x2,P2,[]);
rX = rank(X);

%% Functions

function g = pulse(tau, L, P)
    % Sample Pulse With Delay
    % INPUT:
    % tau : time delay
    % L : dimension
    % P : sampling frequency
    %
    % OUTPUT:
    % g : sampled pulse
    
    % reference pulse
    tri = @(t) triangularPulse(0,2,t);

    % generate g(tau)
    for k=0:(P*L - 1)
            g(k+1) = tri(k*1/P - tau);
    end
end

function h = channel(tau_vec, beta, L, P)
    % Channel Response With Delays
    % INPUT:
    % tau_vec : time delays [tau_1 ... tau_r]
    % beta : 
    % L : dimension
    % P : sampling frequency
    %
    % OUTPUT:
    % h : sampled channel
    
    r = length(beta);
    
    for m = 1:r
        % generate g(t - tau) and h(t)
        tau = tau_vec(m);
        b = beta(m);
        
        g = pulse(tau, L, P);
        
        h_vec(m,:) = b*g;     
    end
    h = sum(h_vec,1);
end



function s = source(N, qpsk_symb)
idxs = randi(length(qpsk_symb),1,N);   % get random index
s(1:N) = qpsk_symb(idxs);
end


function [x, xx] = gen_data(h,s,P,N)

    s_ext = kron(s,[1;zeros(P-1,1)]);
    xx = conv(h,s);

    % sample
    x = xx(1:N*P);
end


