function [xalg,yalg] = runWifiRpca(Network,i,xalg,yalg)

tw = Network.tw;
tr = Network.tr;

if i>tr
     if i>=tw
         tn = tw;
     else 
         tn = i;
     end
     LocCoorTmp = [xalg(:,i-tn+1:i);yalg(:,i-tn+1:i)];
     %[A_hat E_hat iter] = inexact_alm_rpca_stable_fix_rank(LocCoorTmp,1/sqrt(min(size(LocCoorTmp))), 1e-7, 1000,3);
     [A_hat E_hat iter] = inexact_alm_rpca_stable_fix_rank(LocCoorTmp,1,1/sqrt(min(size(LocCoorTmp))), 1e-7, 1000,3);
     xalg(:,i) = A_hat(1:Network.numdev,end);
     yalg(:,i) = A_hat(Network.numdev+1:end,end);
 end