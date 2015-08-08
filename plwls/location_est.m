function [ G_est ] = location_est( G_orient, G_ini, rsstrain, rssloc, loc_set, L)
%Set the search scope and joint localization estmation
%Input "G_orient": the graph construct after orient; "G_ini": graph
%construct by WiFi localization; "rsstrain": rss data in the offline

%"time_second": time second
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


[m, n] = size(G_ini);
G_est = G_orient;
min=inf;

for i=-1/9 : 1/90 : 1/9   %orientation (-20~20 degree)
    G0 = (G_orient - repmat(G_orient(1,:), m, 1))*[cos(i*pi) -sin(i*pi);sin(i*pi) cos(i*pi)];%rotated around the 1st node
    G0 = G0 + repmat(G_orient(1,:), m, 1);
    for j = -2 : 0.1 : 2  %momevent step(-2~2 m, step 0.1m)
        for k = -2 : 0.1 : 2
            A = j*ones(m,1);
            B = k*ones(m,1);
            G = G0 + [A B];
            if max(sum(((G-G_ini).^2),2))<4  %each vertex restricted  inside a circle (set radius as 2m)
                temp = 0;
                idx = knnsearch(loc_set, G, 'k', L);
                for l = 1 : n
                    temp = temp + (sum(rsstrain(idx(l, :), :))/L - rssloc(l, :))*(sum(rsstrain(idx(l, :), :))/L - rssloc(l, :))';
                end
                if temp < min;
                    min = temp;
                    G_est = G;
                end
            end
        end
    end
end
                

end

