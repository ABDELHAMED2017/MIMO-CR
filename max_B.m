function Bmax = max_B(B,H,P,r_min,type)
    [Nrx,Ntx] = size(H);
    lb = exp(r_min)^(1/Ntx);
    C = P/Ntx*eye(Ntx);
    c = norm(C,'fro'); 
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
 x = 1/(a*trace(B*B));
 Bmax = B * x;  
 end
