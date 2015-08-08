function [oneHop,twoHop,oneHopSeed,twoHopSeed]=findNeigh(node,hopDistance,rx)
oneHop1=hopDistance(node,:);
oneHopNode=find(oneHop1==1);
oneHop=oneHopNode(oneHopNode<=rx);
oneHopSeed=oneHopNode(oneHopNode>rx);
twoHopNode=find(oneHop1==2);
twoHop=twoHopNode(twoHopNode<=rx);
twoHopSeed=twoHopNode(twoHopNode>rx);