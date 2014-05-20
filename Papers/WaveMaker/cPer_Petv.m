%%% Compute solitary wave solutions to the Peregrine system
%%% We use a Fourier-type spectral method and Petv. iteration
%%% ----------------------------------------------------------
%%% Author: Denys Dutykh, CNRS - LAMA, University of Savoie
%%% E-mail: Denys.Dutykh@univ-savoie.fr
%%%    URL: http://www.lama.univ-savoie.fr/~dutykh/
%%% ----------------------------------------------------------

close all, clear all
format long e

g = 1.0;                % gravity acceleration
d = 1.0;                % undisturbed water depth

N = 2048;               % number of Fourier modes

l = 40.0;               % half-length of the domain
dx = 2*l/N;             % distance between two points in the real space
x = (1-N/2:N/2)'*dx;    % real space discretization

% vector of wavenumbers
k = [0:N/2 1-N/2:-1]'*pi/l;

j = (N/4+2:N/4*3);      % antialising treatment
k(j) = 0;               % the frequencies we sacrify

c = 1.25;               % solitary wave speed

L = c*(1 + 1/3*(d*k).^2);   % linear operator
Linv = 1./L;                % inverse linear operator

%%% Initial guess specification:
amp = c^2/g - d; k = sqrt(3*amp/(d+amp))/d;
eta0 = amp*sech(0.5*k*x).^2;
u0 = c*eta0./(d+eta0);

fftw('planner', 'patient');

tol = 1e-15;        % iterations stopping criterium
gam = 1.5;          % a parameter in the Petviashvili method
iter = 0;           % iterations count parameter

Err_list = [];

err = inf;
while (err > tol), iter = iter + 1;
    u0_hat = fft(u0);
    eta0 = d*u0./(c-u0);
    Nl = 0.5*u0.^2 + g*eta0;
    Nl_hat = fft(Nl); Nl_hat(j) = 0;
    Lu0_hat = L.*u0_hat;
    nom = abs(sum(u0_hat.*Nl_hat));
    denom = abs(sum(u0_hat.*Lu0_hat));
    S = (denom/nom)^gam;
    u1_hat = S*Linv.*Nl_hat;
    u1 = real(ifft(u1_hat));
    err = norm(u1-u0, inf)
    Err_list = [Err_list; err];
    u0 = u1;
end

eta = d*u0./(c-u0);

max(eta)

figure(1)
subplot(2,1,1)
plot(x, eta, 'b-', 'LineWidth', 1.4), grid on
xlabel('x'); ylabel('\eta(x)')
axis([-l-1 l+1 -0.025 max(eta)+0.05])
title(['Solitary wave solution for c_s = ', num2str(c)])
subplot(2,1,2)
plot(x, u0, 'b-', 'LineWidth', 1.4), grid on
xlabel('x'); ylabel('u(x)')
axis([-l-1 l+1 -0.025 max(u0)+0.05])

sprintf('Iterations number = %d', iter)

figure(2)
semilogy(1:iter, Err_list, 'bo-', 'LineWidth', 1.4), grid on
xlabel('Iteration number'); ylabel('Error');
