%Peer assisted localization algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dist_constr=Network.dis(:,:,20);
rssloc=Network.rssloc(:,:,20);
rsstrain=Network.rsstrain;
G_ini=[X.xcoor(:,20) Y.ycoor(:,20)];
[ G_constr ] = MDS_MAP( dist_constr );
[ G_rotated, dist_constr ] = rotate_longest( G_constr );
[ dire_vec,dire_vec0,dire_vec1,G_rotated0 ] = edge_direct( G_rotated, G_ini,3 );
[ G_orient ] = orientation( dire_vec,dire_vec0,dire_vec1,G_rotated,G_rotated0 );
avg=mean(G_ini-G_orient);
G_orient=G_orient+[avg(1)*ones(length(G_orient(:,1)),1) avg(2)*ones(length(G_orient(:,1)),1)];
[ G_est ] = location_est( G_orient, G_ini, rsstrain, rssloc, 11, 9);

plot(G_ini(:,1),G_ini(:,2),'k+')
hold on
plot(G_rotated(:,1),G_rotated(:,2),'r*')
% plot(G_orient(:,1),G_orient(:,2),'ks')
plot(G_est(:,1),G_est(:,2),'bo')