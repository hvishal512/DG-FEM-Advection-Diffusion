%Evaluate the value of the derivative of a Nth order Lagrange interpolation
%at the intepolation points x_i (GL quad points) for each of the basis
%vectors and return a (N+1)x(N+1) matrix where Dr(n,i)=L_n'(x_i)

% The derivative matrix is later used to calculate the stiffness matrix
function dL = dLagrange(N)
    [Qx,~,~]=lglnodes(N);
    dL = zeros(N+1);
    for n=1:N+1
        dL(n,:) = LagBasisD_n(N,n,Qx);
    end
end

%-----------------------------------------
    
%Evaluate derivative of Lagrange poly basis of point 'n' in an
%interpolation of order N at all of the quadrature points 'Qx'
function [Ld_n] = LagBasisD_n(N,n,Qx)
    NOn=[1:n-1,n+1:N+1]; %The Lagrange bases don't include one for point 'n'
    Ld_n = zeros(N+1,1);
    for i=1:N+1
        NOni=NOn(not(NOn==i)); %The derivative can't include the evaluated point 'i'
        dLagBas_fun = @(x) prod( bsxfun(@minus,x,Qx(NOni)) )./prod( Qx(n)-Qx(NOn) );
        if not(i==n) %Every eval point except the one the basis is built on
            Ld_n(i) = dLagBas_fun(Qx(i));
        elseif i==n %Specially handle the evaluation of the point the basis is built on
            Ld_n(i) = sum( 1./( Qx(i)-Qx(NOn) ) );
        end  
    end
end