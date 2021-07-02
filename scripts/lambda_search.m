close all;
clear all;

% add dependent libs
addpath(genpath('/Users/yuewu/Documents/GitHub/smooth_search/'));
rmpath(genpath('/Users/yuewu/Documents/GitHub/smooth_search/archive'))
addpath(genpath('/Users/yuewu/Documents/GitHub/fda_learn/'));
addpath(genpath('/Users/yuewu/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/matalb.lib/fdaM/'));

% test data simulation 1
rng(1);
x=1:52;
xuse=x/52*2.7;
y_real=xuse.^4-4*xuse.^3+8*xuse+4;
sigma=0.5;
y_raw=y_real+normrnd(0,sigma,[1,length(y_real)]);
y_raw=y_raw.';
x=x.';
y_real=y_real.';
% plot(x,y_raw);
%
methods={'smooth_derivative'};
loglambda_vec= -4:0.25:4;
best_str=struct();
objmat=[];%lambda, true obj, estimated obj with known sigma, GCV
method_vec=[];
for method=methods
  method=method{1};
  bestobj=nan();
  nDer=0;
  % similar to the follow function
  % res=smooth_derivative(y_raw,x,loglambda,nDer,false);
  %
  % smoothing penalty is based on second derivative
  nsmooth=nDer+2;
  % to make use of nth order derivatives, the order of spline basis need at least to be n+2
  norder=nsmooth+2;
  % basis functions = order + #interior knots
  nbasis=length(x)+norder-2;
  bbasis=create_bspline_basis([min(x) max(x)],nbasis,norder,x);
  nlambda=length(loglambda_vec);
  spfdcell=cell(nlambda,1);
  for i=1:nlambda
     lambda=10^loglambda_vec(i);
     fdPari=fdPar(bbasis,int2Lfd(nsmooth),lambda);
     [fdres,dfi,gcvi,~,~,~,y2cMap]=smooth_basis(x,y_raw,fdPari);
     gcv=sum(gcvi(:));
     df=dfi;
     spfdcell{i}=fdres;
     %
     basismat=getbasismatrix(x,getbasis(fdres));
     A_hat=basismat*y2cMap;
     y_smooth=eval_fd(x,fdres);
     obj=obj_sigma(y_raw,full(A_hat),sigma);
     trueobj=sum((y_smooth-y_real).^2)/length(y_raw);
     objmat=[objmat; [lambda trueobj obj gcv]];
     method_vec=[method_vec; method];
  end
end
table(objmat(:,1),objmat(:,2),objmat(:,3),objmat(:,4),method_vec,'VariableNames',{'lambda','true_obj','est_obj','GCV','methods'});
% plot for the best lambda for each mehtod
colorvec={'k','b','r','g','y'};
fig=figure();
plot(x,y_raw,'color',[.6 .6 .6]); title('Function smoothing comparisons');
hold on
for method_i=1:length(methods)
  method=methods{method_i};
  plot(x,best_str.(method).y_smooth,colorvec{method_i},'linewidth',3);
end
legend(['Raw',methods],'location','southeast');
axis tight;
hold off;
