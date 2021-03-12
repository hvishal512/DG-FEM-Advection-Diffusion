% --- Here we attempt to solve the 1D scalar conservation eqn of the form:
%\frac{\partial u)(\partial t} + c\frac{\partial u)(\partial x} -
%alpha\frac{\partial^2 u}{\partial x^2} = 0 with periodic BCs
% u(0,t) = u(1,t) & u_x(0,t) = u_x(1,t)

close all
clear 
tau=2*pi();
c = 0.5;          %Velocity
alpha = 0.05; %Diffusion coefficient

%Rate params
N=3;  % Order of spatial discretization
K=16; % Number of elements

%Discretize the domain into K elements with K+1 nodes, we use a constant
%spacing, but it could be arbitrary
xNode=0:1/K:1;
%@ all the element boundaries, note there are repeats since nodes have
% a coincident at intersections except at the domain boundaries
elemBC = reshape(sort([xNode,xNode(2:end-1)]),2,K)';

%Legendre basis is used since it forms a orthogonal set
[Qx,Qw,~]=lglnodes(N);
%Precompute Legendre values (Vandermonde matrix) for fixed order
%LGL quadrature points
L = zeros(N+1); 
for m=0:N
    temp = legendre(m,Qx)*sqrt((2*m+1)/2); 
    L(:,m+1) = temp(1,:)'; 
end

%--According to Cockburn,Shu 2001 Eqn 2.2 let uh(.,0) be computed by
%\int v \phi = \int u0 \phi for each element (xk-1/2 < x < xk+1/2)
%-We will numerically compute the RHS integral using the quadrature points
BasisWeights = zeros(K,N+1);
map = zeros(K,N+1);
for k=1:K
    % Affine map relating the physical domain [0,1] to comp domain [-1,1]
    map(k,:) = (elemBC(k,2)-elemBC(k,1))*Qx/2 + (elemBC(k,1)+elemBC(k,2))/2;
    for m=0:N
        %Apporimation of RHS integral with gauss quadrature sum
        BasisWeights(k,m+1) = sum(Qw'.*sin(tau*map(k,:)).*L(:,m+1)');
    end
end

% Modal to Nodal conversion for Lagrange interpolation
NodalWeights = L*BasisWeights'; 

%--Now that we have u0 we can begin explicit time stepping with the 
%semi-discrete form of the PDE.
%Flux term for f(u) = cu: g(u-(x),g(u+(x))=u-(x) (Upwind)
%Flux term for f(u) = -alpha*ux: g(u_x-(x),u_x+(x))=u_x+(x) (Downwind)

dT = 0.001; %time step
saveT   = 0.01; %How often do we save the current state for plotting?
endT    = 5;
nsaveT  = floor(saveT/dT);
nT      = floor(endT/dT);
saved   = zeros(N+1,K,(nT/nsaveT)+1);
saved(:,:,1)=NodalWeights;
% Initialize RK4 slopes
k1 = zeros(size(NodalWeights));
k2 = zeros(size(NodalWeights));
k3 = zeros(size(NodalWeights));
k4 = zeros(size(NodalWeights));
j = 1;
for t = 1:1:nT 
    %Iterating through each element 
    for i = 1:K 
        hk = diff(xNode);       %Size of the element
        M = inv(L*L');          %Mass Matrix
        Dr = dLagrange(N)';     %Derivative Matrix
        S = c*((L*L')\Dr) - alpha*(Dr'*Dr);      %Stiffness Matrix
        %S = S - alpha*(Dr*Dr'); %Corrected stiffness matrix due to 
                                %double derivate tem u_{xx}
        flux = zeros(N+1,1);    %Flux vector (upwind)
        dflux = zeros(N+1,1);   %Derivative flux vector (downwind)
        
        u = NodalWeights(:,i);
        source = zeros(size(u));
        for k = 1:size(source,1)
            if u(k) ~= 0
                source(k) = 0.05*exp(-2/abs(u(k)));
            end
        end
        source = source*(1-t/nT);
        source_integral = Qw.*source;
        
        if ~(i == 1 || i == K)  
            ux1 = Dr*NodalWeights(:,i);
            ux2 = Dr*NodalWeights(:,i+1);
            flux(1) = -c*NodalWeights(end,i-1);
            flux(end) = c*NodalWeights(end,i);
            dflux(1) = -alpha*ux1(1);
            dflux(end) = alpha*ux2(1);
            
        %Periodic BCs mean the upwind element for k=1 is k=K,
        elseif i == 1
            ux1 = Dr*NodalWeights(:,i);
            ux2 = Dr*NodalWeights(:,i+1);
            flux(1) = -c*NodalWeights(end,K);
            flux(end) = c*NodalWeights(end,i);
            dflux(1) = -alpha*ux1(1);
            dflux(end) = alpha*ux2(1); 
        %and the downwind element for k=K is k=1
        else
            ux1 = Dr*NodalWeights(:,i);
            ux2 = Dr*NodalWeights(:,1);
            flux(1) = -c*NodalWeights(end,i-1);
            flux(end) = c*NodalWeights(end,i);
            dflux(1) = -alpha*ux1(1);
            dflux(end) = alpha*ux2(1);            
        end
        
        % NET FLUX
        flux = flux - dflux - source_integral ; 
        
        % SEMI-DISCRETE Eqns and RK4 updates
        k1(:,i) =  (2/hk(i))*(M\(S'*NodalWeights(:,i) - flux));
        k2(:,i) =  (2/hk(i))*(M\(S'*(NodalWeights(:,i) + dT*k1(:,i)/2) - flux));
        k3(:,i) =  (2/hk(i))*(M\(S'*(NodalWeights(:,i) + dT*k2(:,i)/2) - flux));
        k4(:,i) =  (2/hk(i))*(M\(S'*(NodalWeights(:,i) + dT*k3(:,i)) - flux));
    end
    %RK4 Update
    NodalWeights = NodalWeights +(dT/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    %Periodically save system state for plotting
    if t/nsaveT==floor(t/nsaveT)    
        j= j+1;
        saved(:,:,j)=NodalWeights;
    end
end

%Plotting the evolution

for i=1:length(saved)
    plot(map',saved(:,:,i))
    axis([0 1 -1.5 1.5])
    text(1.02,1.1,'Time')
    text(1.02,1,num2str((i-1)*saveT));
    title('DG Reactive-Diff Soln','Fontname','Arial','Fontsize',20);
    xlabel('x','Fontname','Arial','Fontsize',16);
    ylabel('u(x)','Fontname','Arial','Fontsize',16);
    pause(0.01);
    %F(i) = getframe(gcf);
end
%video = VideoWriter('with_source_term.avi', 'Uncompressed AVI');
%open(video);
%writeVideo(video,F);
%close(video);
