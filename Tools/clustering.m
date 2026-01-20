function [nmi,ACC,f,Purity,PP]=clustering(S, cls_num, gt)

[C ] = SpectralClustering(S,cls_num);
[nmi] = compute_nmi(gt,C);
ACC = Accuracy(C,double(gt));
[f,~,~] = compute_f(gt,C);
%[~,RI,~,~]=RandIndex(gt,C);




% purity
predLidx = unique(C); 
pred_classnum = length(predLidx);
correnum = 0;
for ci = 1:pred_classnum
    incluster = gt(find(C== predLidx(ci)));
    inclunub = hist(incluster, 1:max(incluster)); 
    if isempty(inclunub) inclunub=0;end;
    correnum = correnum + max(inclunub);
end;
Purity = correnum/length(C);

end
