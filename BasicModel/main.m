
%%function xhat = QuadSys(xhat,uopt,h,n)

% Nonlinear rendezvous docking controller
% of a BlueROV 2 Heavy to a floating WEC (constrained to heave).
% Control of u and w.
% Input u: x and z thrust force
% States eta: [x,0,z,0,0,0]
% States nu: [u,0,w,0,0,0]
% Initial state: x(0) = [0 0 0 0 0 0], u(0) = 0
% Tracking reference: z(t) = sin(0.2*t)
clear 
clc

xhat=zeros(12,1);
uopt = zeros(6,1);

eta = xhat(1:6);
nu = xhat(7:12);
omega = nu(4:6)

tau = zeros(6,1);
tau (1) =  uopt(1);
tau (3) =  uopt(2);

getquadparams

% %state vector: [x,x2] = [eta,nu]
% nu = [u;v;w;p;q;r]; %nu vector
% eta = [x;y;z;phi;theta;psi]; %eta vector


M = mfunc(m,ix,iy,iz); % mass = ridgid body mass + added mass
Minv = inv(M); %mass matrix inverse


h = .01; %time step
l=25; %lenth of simulation seconds
ti = 0 % initial time

for ii =ti+1:h:l+1  %1:n  
    
%     t = (ii-1)*h;
%     S(ii,1) = t;
%     S(ii,2:7) = transpose(eta); % eta'
%     S(ii,8:13) = transpose(nu);

            %call matrix building functions
            G = gfunc(m); %hydrostatic forces
            J = jfunc(eta(4),eta(5),eta(6)); %refernce frame transform matrix
            C = cfunc(omega, M);
            
            xhat = eulerfunc2(h,nu,eta,Minv,G,J,tau,M,C);

end


% end