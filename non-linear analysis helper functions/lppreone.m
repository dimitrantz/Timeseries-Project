function [xpre,neiindV] = lppreone(xV,TreeRoot,winnowV,m,tau,nnei,q,mindist,tarintree)
% xpre = lppreone(xV,TreeRoot,winnowV,m,tau,nnei,q,mindist,tarintree)
% 
f = 1.2;  % factor to increase the distance if not enough neighbors are found 
tarV = winnowV((m-1)*tau+1:-tau:1)';
distnow = mindist;
if tarintree
    [neiM,neidisV,neiindV]=kdrangequery(TreeRoot,tarV,distnow);
    while length(neiindV)<nnei+1
        distnow = f*distnow;
        [neiM,neidisV,neiindV]=kdrangequery(TreeRoot,tarV,distnow);
    end
    [oneidisV,oneiindV]=sort(neidisV);
    neiindV = neiindV(oneiindV(2:nnei+1));
    neidisV = neidisV(oneiindV(2:nnei+1));
    neiM = neiM(oneiindV(2:nnei+1),:);
    yV = xV(neiindV+(m-1)*tau+1);
    if q==0 | nnei==1
        xpre = mean(yV);
    else
        mneiV = mean(neiM);
        my = mean(yV);
        zM = neiM - ones(nnei,1)*mneiV;
        [Ux, Sx, Vx] = svd(zM, 0);
        tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
        lsbV = tmpM * (yV - my);
        xpre = my + (tarV-mneiV) * lsbV;
    end
elseif nnei==1
    [neiindV,neidisV,TreeRoot]=kdtreeidx([],tarV,TreeRoot);
    xpre = xV(neiindV+(m-1)*tau+1);
else
    [neiM,neidisV,neiindV]=kdrangequery(TreeRoot,tarV,distnow);
    while length(neiindV)<nnei
        distnow = f*distnow;
        [neiM,neidisV,neiindV]=kdrangequery(TreeRoot,tarV,distnow);
    end
    [oneidisV,oneiindV]=sort(neidisV);
    neiindV = neiindV(oneiindV(1:nnei));
    neidisV = neidisV(oneiindV(1:nnei));
    neiM = neiM(oneiindV(1:nnei),:);
    yV = xV(neiindV+(m-1)*tau+1);
    if q==0 
        xpre = mean(yV);
    else
        mneiV = mean(neiM);
        my = mean(yV);
        zM = neiM - ones(nnei,1)*mneiV;
        [Ux, Sx, Vx] = svd(zM, 0);
        tmpM = Vx(:,1:q) * inv(Sx(1:q,1:q)) * Ux(:,1:q)';
        lsbV = tmpM * (yV - my);
        xpre = my + (tarV-mneiV) * lsbV;
    end
end