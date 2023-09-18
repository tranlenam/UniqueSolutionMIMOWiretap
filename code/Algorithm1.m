function [SecCapAPGM,obj_seq] = Algorithm1(Hb,He,nA,nB,nE,maxIter,P0)
gammau=2;
obj_seq =zeros(maxIter,1);
X_past=  ones(nA,nA);%x0
Y_now = X_past;%y1
t=0.5;
myalpha=0.2;
L0 = 0.01;
mybeta=L0;
for iTer=1:maxIter
    gradObj=(Hb'*(inv(eye(nB)+Hb*Y_now*Hb'))*Hb)-He'*(inv(eye(nE)+He*Y_now*He'))*He;

    while(1)
        X_temp=projectSPC(Y_now+1/mybeta*gradObj,P0);
        F = ComputeSecrecyRate(Hb,He,X_temp);
        Q = ComputeSecrecyRate(Hb,He,Y_now)...
            +real(trace((X_temp-Y_now)*gradObj))...
            -mybeta/2*norm(X_temp-Y_now,'fro')^2;
        if(F<Q)
            mybeta=mybeta*gammau;
        else
            break
        end
    end
    mybeta = max(L0,mybeta/gammau);

    X_now=X_temp;

    V_now=X_now+myalpha*(X_now-X_past);

    F_X=ComputeSecrecyRate(Hb,He,X_now) ;
    F_v= ComputeSecrecyRate(Hb,He,V_now) ;

    [flag1] = IsFeasible(V_now,P0);
    if (F_v>F_X)&&(flag1==1)
        Y_next=V_now;
        myalpha=min(myalpha/t,1);
    else
        Y_next=X_now;
        myalpha=t*myalpha;
    end

    %
    X_past=X_now;
    Y_now=Y_next;

    obj_seq(iTer) = ComputeSecrecyRate(Hb,He,X_now); 
    if (iTer>5)
        if(abs(obj_seq(iTer)-obj_seq(iTer-5))<=1e-5)
            break
        end
    end

end
obj_seq(iTer+1:end)=[];
SecCapAPGM=obj_seq(end);

end

