function out=calerr_minflag(D,EstD,minFlag)
D=max(D,0);
stress=sqrt((sum(sum((D-EstD).^2)))/(sum(sum(D.^2))));
% re=(abs(D-EstD))./(D+0.1);
switch minFlag
    case 'min'
    re=(abs(D-EstD))./min(D,EstD);
    case 'nomin'
    re=(abs(D-EstD))./D;
    otherwise
        disp('please choose min or nomin!');
end
% save(['re.mat'],'re');
% for i=1:N
%     re(i,i)=0;
% end
% figure;
% pcolor(re);
D1=D(:);
tmp=find(D<=0);
re1=re(:);
re1(tmp)=[];
E=sort(re1);
mre=mean(E);
N1=length(E);
%max(E)
npre=E(ceil(0.9*N1)); %0.9为npre
median=E(ceil(0.5*N1));%0.5为中值
total_re=[];
output_re = store_re(re1', 1, 1000);
total_re=[total_re;output_re];

out.npre = npre;
out.stress = stress;
out.median = median;
out.re1 = re1;
out.mre = mre;

% figure;
% h1 = plot(0:1/1000:1-1/1000, output_re, 'b--');
% set(h1, 'LineWidth', 2);hold on;
% h0 = plot(0:1, [0.9 0.9], 'r:');hold on;
% xlabel('Relative Error', 'FontSize', 16);
% ylabel('Cumulative Distribution Function', 'FontSize', 16);axis([0 1 0 1]);