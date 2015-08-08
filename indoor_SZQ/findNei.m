function tNei=findNei(node,hopDistance,N)
Nei=hopDistance(node,:);
tNei=find(Nei<=N);
