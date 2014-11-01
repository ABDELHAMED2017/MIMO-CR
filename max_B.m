function Bmax = max_B(B,H,P,r_min,type)
    [Nrx,Ntx] = size(H);
    lb = exp(r_min)^(1/Ntx);
    C = P/Ntx*eye(Ntx);
    c = norm(C,'fro'); 
    B = B/ norm(B,'fro');
    switch type
        case 'sum-power'   
            cvx_begin sdp quiet
                variable Si(Nrx,Nrx) symmetric
                variable a
                minimize ( a )
                subject to  
                det_rootn(eye(Ntx) + H'*Si*H) >= lb
                0.5*(Si + Si') >= 0
                0.5*(Si + Si') <= 0.5/P*Ntx*c*c*a*(B + B')
            cvx_end
%         B1 = eye(Nrx);
%         cvx_begin sdp
%             variable Si(Nrx,Nrx) symmetric
%             variable Y1(Nrx,Nrx) symmetric
%             variable Y2(Nrx,Nrx) symmetric
%             variable a
%             variable b
%             variable la
%             minimize ( a )
%             subject to 
%             det_rootn(eye(Ntx) + H'*Si*H) >= lb
%             0.5*(Si + Si') >= 0
%             0.5*(Si + Si') <= 0.5/P*Ntx*c*c*a*(B1+B1') + Y1
%             trace(B1*Y1) + trace(B*Y2) == 0
%             0.5*(Si + Si') <= 0.5/P*Ntx*c*c*a*(B + B')+ Y2  
%             Y2 == la*eye(Nrx)
%         cvx_end      
%         
%         b = 1/a - trace(B1*B1);
%         Bmax = sqrt(b)*B
%         s = 1/(trace(B1*B1) + trace(Bmax*Bmax));
%         cvx_begin sdp
%             variable Si(Nrx,Nrx) symmetric
%             variable Y1(Nrx,Nrx) symmetric
%             variable Y2(Nrx,Nrx) symmetric
%             variable la
%             maximize ( det_rootn(eye(Ntx) + H'*Si*H) )
%             subject to 
%             0.5*(Si + Si') >= 0
%             0.5*(Si + Si') <= 0.5/P*Ntx*c*c*s*(B1+B1') + Y1
%             trace(B1*Y1) + trace(Bmax*Y2) == 0
%             0.5*(Si + Si') <= 0.5/P*Ntx*c*c*s*(Bmax + Bmax')+ Y2  
%             Y2 == la*eye(Nrx)
%         cvx_end
%         normbsq = 1/a - trace(B1*B1);
%         Bmax = B * sqrt(normbsq)
%         s = trace(C*C) / (trace(B1*B1) + trace(Bmax*Bmax))
%         cvx_begin sdp
%             variable Si(Nrx,Nrx) symmetric
%             variable Y1(Nrx,Nrx) symmetric
%             variable Y2(Nrx,Nrx) symmetric
%             variable la
%             dual variable R1
%             dual variable R2
%             maximize ( det_rootn(eye(Ntx) + H'*Si*H) )
%             subject to  
%             0.5*(Si + Si') >= 0
%             R1: 0.5*(Si + Si') <= 0.5/P*Ntx*s*(B1+B1') + Y1
%             trace(B1*Y1) + trace(B*Y2) == 0
%             R2: 0.5*(Si + Si') <= 0.5/P*Ntx*s*(B + B') + Y2  
%             Y2 == la*eye(Nrx)
%         cvx_end
%         log(det(eye(Ntx) + H'*Si*H))

        case 'per-antenna'
            cvx_begin sdp quiet
                variable Si(Nrx,Nrx) symmetric
                variable Y(Nrx,Nrx) symmetric
                variable a
                minimize ( a )
                subject to  
                det_rootn(eye(Ntx) + H'*Si*H) >= lb
                0.5*(Si + Si') >= 0
                0.5*(Si + Si') <= 0.5/P*Ntx*c*c*a*(B + B') + Y
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
                variable a
                minimize ( a )
                subject to  
                det_rootn(eye(Ntx) + H'*Si*H) >= lb
                0.5*(Si + Si') >= 0
                0.5*(Si + Si') <= 0.5/P*Ntx*c*c*a*(B + B') + Y
                trace(B*Y) == 0
            cvx_end  
    end    
 b = 1/a;
 Bmax = B * b;  
 end