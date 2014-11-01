function [rate,R,Q] = min_max_dl(H,C,B,type)
[Nrx,Ntx] = size(H);
P = trace(C);
% Compute the worst case noise
switch type
    case 'sum-power'    
          cvx_begin sdp quiet
            variable L(Ntx,Ntx)
            variable A(Nrx,Nrx)
            variable Yp(Nrx,Nrx) symmetric
            variable muh
            minimize ( -det_rootn(L) )
            subject to 
            trace(L) - Ntx + muh*P == 0
            A >= 0.5*(H*L*H' + (H*L*H')')
            muh >= 0
            0.5*(L + L') >= 0
            0.5*(A + A') >= 0           
            A <= muh*B + Yp
            trace(B*Yp) == 0
        cvx_end
%         cvx_begin sdp quiet
%             variable L(Ntx,Ntx)
%             variable A(Nrx,Nrx) symmetric
%             variable A1(Nrx,Nrx) symmetric
%             variable A2(Nrx,Nrx) symmetric
%             variable Yp(Nrx,Nrx) symmetric
%             variable muh
%             minimize ( -det_rootn(L) )
%             subject to 
%             trace(L) - Ntx + muh*P == 0
%             A == A1 + A2
%             A >= 0.5*(H*L*H' + (H*L*H')')
%             muh >= 0
%             0.5*(L + L') >= 0
%             0.5*(A + A') >= 0
%             0.5*(A1 + A1') >= 0
%             0.5*(A2 + A2') >= 0
%             A2 == muh*B + Yp
%             trace(B*Yp) == 0
%             A1 == muh*eye(Nrx)
%         cvx_end
        R = A/muh;
    case 'per-antenna'
        cvx_begin sdp quiet
            variable L(Ntx,Ntx)
            variable A(Nrx,Nrx)
            variable Yp(Nrx,Nrx) symmetric
            variable muh
            minimize ( -det_rootn(L) )
            subject to 
            trace(L) - Ntx + muh*P == 0
            A >= 0.5*(H*L*H' + (H*L*H')')
            muh >= 0
            0.5*(L + L') >= 0
            0.5*(A + A') >= 0           
            A <= muh*B + Yp
            for n=1:Nrx
                Yp(n,n) == 0
            end
        cvx_end
        R = A/muh;
    case 'shape'
        R = B;
end       
    
Hd = R^-0.5*H;
% compute the transmit covariance
cvx_begin sdp quiet
    variable Q(Ntx,Ntx) symmetric
    variable Z(Ntx,Ntx) symmetric
    maximize ( det_rootn(eye(Nrx) + Hd*Q*Hd')  )
    subject to  
    Q <= C + Z
    trace(C*Z) == 0
    Q >= 0
cvx_end
rate = log(det(eye(Nrx) + Hd*Q*Hd'));