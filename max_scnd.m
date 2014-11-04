function [rate,Q,A] = max_scnd(H,Hi,P,B,type)
[Nrx,Ntx] = size(H);
Irx = size(Hi,1);

switch type
    case 'sum-power'   
        cvx_begin sdp quiet
            variable Q(Ntx,Ntx) symmetric
            variable A(Irx,Irx) symmetric
            variable Y(Irx,Irx) symmetric
            maximize ( det_rootn(eye(Nrx) + H*Q*H') )
            subject to  
            0.5*(Q + Q') >= 0
            trace(Q) <= P
            A == eye(Irx) + Hi*Q*Hi' 
            0.5*(A + A') <= 0.5*(B + B') + Y
            trace(B*Y) == 0
        cvx_end
       
    case 'per-antenna'
        cvx_begin sdp quiet
            variable Q(Ntx,Ntx) symmetric
            variable A(Irx,Irx) symmetric
            variable Y(Irx,Irx) symmetric
            maximize ( det_rootn(eye(Nrx) + H*Q*H') )
            subject to  
            0.5*(Q + Q') >= 0
            trace(Q) <= P
            A == eye(Irx) + Hi*Q*Hi' 
            0.5*(A + A') <= 0.5*(B + B') + Y
            for n = 1:Irx
                Y(n,n) == 0
            end
        cvx_end
        
    case 'shape'
        cvx_begin sdp quiet
            variable Q(Ntx,Ntx) symmetric
            variable A(Irx,Irx) symmetric
            maximize ( det_rootn(eye(Nrx) + H*Q*H') )
            subject to  
            0.5*(Q + Q') >= 0
            trace(Q) <= P
            A == eye(Irx) + Hi*Q*Hi' 
            0.5*(A + A') <= 0.5*(B + B') 
        cvx_end
        
end       
rate = log(det(eye(Nrx) + H*Q*H'));