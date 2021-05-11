%% SP for COMM HW 1 Q 2
% Convolutive Model

%% Question 1 (a)
L = 3;
syms t
t = 0:.1:L;

tres = .1;                          % time delay resolution
tau_array = (0:tres:1-tres);        % delay values
idx=randperm(length(tau_array),1);  % get random index
tau = tau_array(idx);               % get random delay 'tau'

% Alternatively, choose a specific delay value by commenting the tau above
% and entering a value below:
% tau = 0.5;    % choose from 0, 0.1, ... 0.8, 0.9

% choose arbitrary P, L values for now
P = 1.5*(1/tres);
L = 3;

% sample pulse with delay
g = pulse(t,tau,L,P);

% new time axis
t22 = linspace(0,3,L*P);
% % plot over-sampled g(tau)
% figure(1); plot(t2,g);
% title('sampled pulse'); xlabel('t');
% 
% % plot reference pulse and delayed pulse
% for k = 1
%     tri = @(t) triangularPulse(0,2,t);
%     figure(1);
%     plot(t,tri(t));
%     hold on;
%     plot(t22,g);
%     ylabel('g(t)'); xlabel('t');
%     str = {'Time Delay tau',num2str(tau)};
%     annotation('textbox',[0.2 0.5 0.3 0.3],'String',str,'FitBoxToText','on');
%     title('reference pulse and delayed pulse');
% end

%% Question 1(b)
r = 2;                                  % arbitrary number of paths
idxs = randperm(length(tau_array),r);   % get random index
tau_vec = tau_array(idxs);              % generate time delay vector

beta = 1;                               % gains, column vector

h = channel(t,tau_vec, beta, L, P);

%% Question 1 (c)
qpsk_symb = [1 1j -1 -1j];
N = 8;
s = source(N, qpsk_symb)';

%% Question 1 (d)
x = gen_data(h,s,P,N);

%% Question 2 - PROBLEM

r2 = 2;
tau2 = [0.1 0.6]';
phi1 = 0;   % deg
phi2 = 15;  % deg
beta2 = [1*exp(1j*deg2rad(phi1)) 0.7*exp(1j*deg2rad(phi2))]';
P2 = 5;
N2 = 50;

t2 = linspace(1,N2,N2);
L2 = N2+2;

h2 = channel(t2,tau2, beta2, L2, P2)';

s2 = source(N2, qpsk_symb)';

[x2, xx2] = gen_data(h2,s2,P2,N2);

% plot results
lenx = length(x2);
t_ax2 = linspace(0,lenx+1,lenx);
figure(2); plot(t_ax2, real(x2)); hold on; plot(t_ax2, imag(x2))
legend('real','imag');title('real and imaginary parts of x2')
xlim([0 300]);

lenh = length(h2);
t_ax3 = linspace(0,lenh+1,lenh);
figure(3); plot(t_ax3, real(h2)); hold on; plot(t_ax3, imag(h2))
legend('real','imag');title('real and imaginary parts of h2')
xlim([0 15]);

%% Question 3
X = reshape(x2,P2,[]);

%% Functions

function g = pulse(t,tau, L, P)
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
    
    % Alternative
    %g = resample(tri(t-tau),P,1);
end

function h = channel(t,tau_vec, beta, L, P)
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
    T = length(t);
    
    for m = 1:r
        % generate g(t - tau) and h(t)
        tau = tau_vec(m);
        b = beta(m);
        
        g = pulse(t,tau, L, P);
        
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
    x = resample(xx,P,1);
    %x = xx(1:1/P:(P*N - 1)*1/P);
end


