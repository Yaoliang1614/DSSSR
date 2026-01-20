clear;
clc;
load 'bbcsport_2view.mat';
AA=[];
v=size(data,2);
for i = 1:v
    data{i}=full(data{i}');
    data{i} = NormalizeFea(data{i},0);
    AA=[AA;data{i}];
end
gnd=truth;
Yfea=AA;

%% Subspace segmentation
Par.rho = 1.9;
Par.maxIter=50;
nCluster = length( unique( gnd ) ) ;
Par.lambda_1 = 2;
Par.lambda_2 = 10;
Par.s = 0.1;
for kk=1:15
    [ C,Z]= DS_SSRLSR( Yfea , Par ) ;
    C= (  C +  C'  )/2;
    [NMI_c(kk),ACC_c(kk),F_c(kk)]=clustering(C, nCluster, gnd);
    fprintf('iter=%dth  C: NMI = %f, ACC = %f, F = %f \n',kk,NMI_c(kk),ACC_c(kk),F_c(kk));
end
ACC=mean(acc);
NMI=mean(nmi);
F=mean(f);
Purity=mean(purity);
[ACC,NMI,F,Purity]