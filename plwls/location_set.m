function [ loc_set ] = location_set( L, L0 )
%location train, where L, L0 represent the size of train area 
%return the loc coor of train location

loc_set = zeros(L * L0, 2);
k = 1;
for i = 1 : L
    for j = 1 : L0
        loc_set(k, :) = [i-1 j-1];
        k = k + 1;
    end
end

end

