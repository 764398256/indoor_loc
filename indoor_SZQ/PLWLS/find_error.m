function [ idx ] = find_error( cu_num, n, percent )
% ����ٷֱ����ֵ
%   �˴���ʾ��ϸ˵��

L = length(cu_num);
for i = 1:L
    if(sum(cu_num(1:i))>n*percent)
        break;
    end
end

idx = i-1;


end

