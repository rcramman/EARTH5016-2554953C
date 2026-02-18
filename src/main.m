%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************


%*****  Initialise Model Setup

% create coordinate vectors
xc = h/2:h:W-h/2;    % coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;
xf = 0:h:W;            % coordinate vectore for cell face positions [m]
zf = 0:h:D;
[Xc,Zc] = meshgrid(xc,zc);


% set up ghosted index lists for boundary conditions
switch BC
    case 'periodic'

        % example periodic indexing for N=4
        % 3-point stencil |-- i-1 --|-- i --|-- i+1 --|
        % 5-point stencil |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
        ix3 = [     Nx,1:Nx,1 ];
        ix5 = [Nx-1,Nx,1:Nx,1,2];
        iz3 = [     Nz,1:Nz,1 ];
        iz5 = [Nz-1,Nz,1:Nz,1,2];
    case 'insulating'
        % example non-periodic indexing for N=4 
        % 3-point stencil            |-- i-1 --|-- i --|-- i+1 --|
        % 5-point stencil  |-- i-2 --|-- i-1 --|-- i --|-- i+1 --|-- i+2 --|
        ix3 = [   1,1:Nx,Nx     ];
        ix5 = [ 1,1,1:Nx,Nx,Nx  ];  
        iz3 = [   1,1:Nz,Nz     ];
        iz5 = [ 1,1,1:Nz,Nz,Nz  ];
end

% set initial coefficient fields
kT = kT0 .* ones(Nz,Nx);
cP = cP0 .* ones(Nz,Nx);
rho = rho0 .* ones(Nz,Nx);
Qr = Qr0 .* ones(Nz,Nx);

% set initial velocity field
w = w0 .* ones(Nz+1,Nx);
u = u0 .* ones(Nz,Nx+1);

% Initialise time count variables
t = 0;  % initial time [s]
k = 0;  % initial time step count

% set time step size
dt_adv = (h/2) / (max(abs(u(:))) + max(abs(w(:)))+eps);
dt_dff = (h/2)^2/(max(kT(:)./rho(:)./cP(:))+eps);
dt     = CFL * min(dt_adv,dt_dff); % time step [s]

% set initial condition for temperature at cell centres
sgmt = sqrt(sgm0^2 + 2*k0*t);
T   = T0 + dT .* exp(-((Xc-W/2  -u0*t).^2 + (Zc-D/2-w0*t  ).^2)./(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2-u0*t  ).^2 + (Zc-D/2-D-w0*t).^2)./(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2-u0*t  ).^2 + (Zc-D/2+D-w0*t).^2.)/(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2-W-u0*t).^2 + (Zc-D/2-w0*t  ).^2.)/(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2-W-u0*t).^2 + (Zc-D/2-D-w0*t).^2.)/(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2-W-u0*t).^2 + (Zc-D/2+D-w0*t).^2.)/(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2+W-u0*t).^2 + (Zc-D/2-w0*t  ).^2.)/(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2+W-u0*t).^2 + (Zc-D/2-D-w0*t).^2.)/(4*sgmt^2)) ...
         + dT .* exp(-((Xc-W/2+W-u0*t).^2 + (Zc-D/2+D-w0*t).^2.)/(4*sgmt^2));
Tin = T;                                         % store initial condition for plotting
Ta  = T; 

% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,Ta,0); %Tin


%*****  Solve Model Equations

while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

    % calculate T-dependent density
    %rho = rho0.*(1-a_T*(T-T0));

    % select time integration scheme
    switch TINT
        case 'FE1'  % 1st-order Forward Euler time integration scheme
            
            % get rate of change
            dTdt = advection(T,h,w,u,ix5,iz5,ADVN) ...
                   + diffusion(T,kT,h,ix3,iz3)./(rho.*cP) ...
                   + Qr./(rho.*cP);

        case 'RK2'  % 2nd-order Runge-Kutta time integration scheme
            
            dTdt_half = diffusion(T,kT,h,ix3,iz3)./(rho.*cP) ...
                      + advection(T,h,w,u,ix5,iz5,ADVN) ...
                      + Qr./(rho.*cP);
            dTdt      = diffusion(T+dTdt_half.*dt/2,kT,h,ix3,iz3)./(rho.*cP) ...
                      + advection(T+dTdt_half.*dt/2,h,w,u,ix5,iz5,ADVN) ...
                      + Qr./(rho.*cP);

    end

    % update temperature
    T = T + dTdt * dt;

    % get analytical solution at time t
    sgmt = sqrt(sgm0^2 + 2*k*t);
    Ta   = T0 + dT .* exp(-((Xc-W/2  -u0*t).^2 + (Zc-D/2-w0*t  ).^2)./(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2-u0*t  ).^2 + (Zc-D/2-D-w0*t).^2)./(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2-u0*t  ).^2 + (Zc-D/2+D-w0*t).^2.)/(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2-W-u0*t).^2 + (Zc-D/2-w0*t  ).^2.)/(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2-W-u0*t).^2 + (Zc-D/2-D-w0*t).^2.)/(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2-W-u0*t).^2 + (Zc-D/2+D-w0*t).^2.)/(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2+W-u0*t).^2 + (Zc-D/2-w0*t  ).^2.)/(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2+W-u0*t).^2 + (Zc-D/2-D-w0*t).^2.)/(4*sgmt^2)) ...
              + dT .* exp(-((Xc-W/2+W-u0*t).^2 + (Zc-D/2+D-w0*t).^2.)/(4*sgmt^2));



    % plot model progress
    if ~mod(k,nop)
        makefig(xc,T,Tin,Ta,t/yr);
        pause(0.1);
    end

end


%*****  calculate and display numerical error norm
Err = rms(T-Ta)/rms((Ta));

disp(' ');
disp(['Advection scheme: ',ADVN]);
disp(['Time integration scheme: ',TINT]);
disp(['Numerical error = ',num2str(Err)]);
disp(' ');



%*****  Utility Functions  ************************************************

%*****  Function to make output figure

function makefig(x,z,T,Ta,t) %Tin

subplot(2,1,1)
imagesc(x,z,T); axis equal tight; colorbar

ylabel('z [m]','FontSize',15)
title(['Temperature [C]; time = ',num2str(t),' [yr]'],'FontSize',17)

subplot(2,1,2)
imagesc(x,z,T-Ta); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Num. Error [C]','FontSize',17)

drawnow;

end


%*****  Function to calculate diffusion rate

function dfdt = diffusion(f,k,h,ix,iz)

% input arguments
% f:    diffusing scalar field
% k:    diffusion coefficient
% dx:   grid spacing
% ind:  ghosted index list

% output variables
% dfdt: diffusion rate of scalar field f

% calculate diffusive flux coefficient at cell faces
kfz = (k(iz(1:end-1),:)+k(iz(2:end),:))/2;
kfx = (k(:,ix(1:end-1))+k(:,ix(2:end)))/2;

% calculate diffusive flux of scalar field f
qz = - kfz .* diff(f(iz,:),1,1)/h;
qx = - kfx .* diff(f(:,ix),1,2)/h;

% calculate diffusion flux balance for rate of change
dfdt = - diff(qz,1,1)/h ...
       - diff(qx,1,2)/h;


end


%*****  Function to calculate advection rate

function dfdt = advection(f,h,w,u,ix,iz,ADVN)

% input arguments
% f:    advected scalar field
% u:    advection velocity
% dx:   grid spacing
% ind:  ghosted index list
% ADVN: advection scheme

% output variables
% dfdt: advection rate of scalar field

% split the velocities into positive and negative
u_pos = max(0,u);    % positive velocity (to the right)
u_neg = min(0,u);    % negative velocity (to the left)
w_pos = max(0,w);
w_neg = min(0,w);

% get values on stencil nodes
f_imm  = f(:,ix(1:end-4));  % i-2
f_im   = f(:,ix(2:end-3));  % i-1
f_ic   = f(:,ix(3:end-2));  % i
f_ip   = f(:,ix(4:end-1));  % i+1
f_ipp  = f(:,ix(5:end  ));  % i+2

f_jmm  = f(iz(1:end-4),:);  % j-2
f_jm   = f(iz(2:end-3),:);  % j-1
f_jc   = f(iz(3:end-2),:);  % j
f_jp   = f(iz(4:end-1),:);  % j+1
f_jpp  = f(iz(5:end  ),:);  % j+2

% get interpolated field values on i+1/2, i-1/2 cell faces
switch ADVN
    case 'UPW1'   % 1st-order upwind scheme
        % positive velocity
        f_ip_pos = f_ic;     % i+1/2
        f_im_pos = f_im;     % i-1/2
        f_jp_pos = f_jc;
        f_jm_pos = f_jm;

        % negative velocity
        f_ip_neg = f_ip;     % i+1/2
        f_im_neg = f_ic;     % i-1/2
        f_jp_neg = f_jp;
        f_jm_neg = f_jc;

    case 'CFD2'  % 2nd-order centred finite-difference scheme
        % positive velocity
        f_ip_pos = (f_ic+f_ip)./2;     % i+1/2
        f_im_pos = (f_ic+f_im)./2;     % i-1/2
        f_jp_pos = (f_jc+f_jp)./2;
        f_jm_pos = (f_jc+f_jm)./2;

        % negative velocity
        f_ip_neg = f_ip_pos;          % i+1/2
        f_im_neg = f_im_pos;          % i-1/2
        f_jp_neg = f_jp_pos;         
        f_jm_neg = f_jm_pos;          

    case 'UPW3'  % 3rd-order upwind scheme
        % positive velocity
        f_ip_pos = (2*f_ip + 5*f_ic - f_im )./6;     % i+1/2
        f_im_pos = (2*f_ic + 5*f_im - f_imm)./6;     % i-1/2

        f_jp_pos = (2*f_jp + 5*f_jc - f_jm )./6;     % i+1/2
        f_jm_pos = (2*f_jc + 5*f_jm - f_jmm)./6;     % i-1/2

        % negative velocity
        f_ip_neg = (2*f_ic + 5*f_ip - f_ipp)./6;     % i+1/2
        f_im_neg = (2*f_im + 5*f_ic - f_ip )./6;     % i-1/2

        f_jp_neg = (2*f_jc + 5*f_jp - f_jpp)./6;     % i+1/2
        f_jm_neg = (2*f_jm + 5*f_jc - f_jp )./6;     % i-1/2
end

% calculate advection fluxes on i+1/2, i-1/2 cell faces

% positive velocity
q_ip_pos = u_pos(:,2:end  ).*f_ip_pos;
q_im_pos = u_pos(:,1:end-1).*f_im_pos;

q_jp_pos = w_pos(2:end,:  ).*f_jp_pos;
q_jm_pos = w_pos(1:end-1,:).*f_jm_pos;

% negative velocity
q_ip_neg = u_neg(:,2:end).*f_ip_neg;
q_im_neg = u_neg(:,1:end-1).*f_im_neg;

q_jp_neg = w_neg(2:end,:).*f_jp_neg;
q_jm_neg = w_neg(1:end-1,:).*f_jm_neg;

% advection flux balance for rate of change
div_q_pos = (q_ip_pos - q_im_pos)/h + (q_jp_pos - q_jm_pos)/h;  % positive velocity
div_q_neg = (q_ip_neg - q_im_neg)/h + (q_jp_neg - q_jm_neg)/h;  % negative velocity

div_q     = div_q_pos + div_q_neg;     % combined
dfdt      = - div_q;                   % advection rate

end


%*****  Function to calculate darcy flux

%function [w,u] = darcyflux(k,h,rho,ix,iz,g)

%qz = - k .* diff(f(),1,1)/h;
%qx = - k .* diff(f(:,ix),1,2)/h;






