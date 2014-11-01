function [rate,Om,Si] = min_max_ul(H,C,B,type)
[Nrx,Ntx] = size(H);
Hu = H*C^-0.5;
c = norm(C,'fro');
b = norm(B,'fro');


B1 = eye(Nrx);
B2 = B;
s = trace(C*C)/(trace(B1*B1) + trace(B2*B2));
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
%         cvx_begin sdp
%             variable Si(Nrx,Nrx) symmetric
%             variable Y1(Nrx,Nrx) symmetric
%             variable Y2(Nrx,Nrx) symmetric
%             variable la
%             dual variable R1
%             dual variable R2
%             maximize ( det_rootn(eye(Ntx) + Hu'*Si*Hu) )
%             subject to  
%             0.5*(Si + Si') >= 0
%             R1: 0.5*(Si + Si') <= 0.5*s*(B1+B1') + Y1
%             trace(B1*Y1) + trace(B2*Y2) == 0
%             R2: 0.5*(Si + Si') <= 0.5*s*(B2 + B2') + Y2  
%             Y2 == la*eye(Nrx)
%         cvx_end
%         cvx_begin sdp
%             variable Si(Nrx,Nrx) symmetric
%             variable Y1(Nrx,Nrx) symmetric
%             variable Y2(Nrx,Nrx) symmetric
%             variable la
%             dual variable R1
%             dual variable R2
%             maximize ( det_rootn(eye(Ntx) + Hu'*Si*Hu) )
%             subject to  
%             0.5*(Si + Si') >= 0
%             R1: 0.5*(Si + Si') <= 0.5*s*(B1+B1') + Y1
%             trace(B1*Y1) + trace(B2*Y2) == 0
%             R2: 0.5*(Si + Si') <= 0.5*s*(B2 + B2') + Y2  
%             Y2 == la*eye(Nrx)
%         cvx_end

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