function [rate,Om,Si] = min_max_ul(H,C,B,type)
[Nrx,Ntx] = size(H);
Hu = H*C^-0.5;
c = norm(C,'fro');
b = norm(B,'fro');



s = trace(C*C)/trace(B*B) ;
% Compute the worst case noise
switch type
    case 'sum-power'   
        cvx_begin sdp quiet
            variable Si(Nrx,Nrx) symmetric
            maximize ( det_rootn(eye(Ntx) + Hu'*Si*Hu) )
            subject to  
            0.5*(Si + Si') >= 0
            0.5*(Si + Si') <= 0.5*c*c/b/b*(B + B') 
        cvx_end

    case 'per-antenna'
        cvx_begin sdp quiet
            variable Si(Nrx,Nrx) symmetric
            variable Y(Nrx,Nrx) symmetric
            maximize ( det_rootn(eye(Ntx) + Hu'*Si*Hu) )
            subject to  
            0.5*(Si + Si') >= 0
            0.5*(Si + Si') <= 0.5*c*c/b/b*(B + B') + Y
            for n = 1:Nrx
                for m = 1:Nrx
                    if n ~= m
                        Y(n,m) == 0
                    end
                end
            end
            trace(B*Y) == 0
        cvx_end
        
    case 'shape'
        cvx_begin sdp quiet
            variable Si(Nrx,Nrx) symmetric
            variable Y(Nrx,Nrx) symmetric
            maximize ( det_rootn(eye(Ntx) + Hu'*Si*Hu) )
            subject to  
            0.5*(Si + Si') >= 0
            0.5*(Si + Si') <= 0.5*c*c/b/b*(B + B') + Y
            trace(B*Y) == 0
        cvx_end
        
end       
Om = C;
Si = C^-0.5*Si*C^-0.5;
rate = log(det(eye(Nrx) + H'*Si*H));
