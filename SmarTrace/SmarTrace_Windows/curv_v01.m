function curv_v01(filename)

%reading sampled chains data from SmTr1.m
dat = load(filename);
FN = strrep(filename(1:end-4),'_','-');
mkdir(FN);
data = dat.sampled_segments;
cr = [];
chain_number = [];
l_on_chain = [];

for i=1:length(data)
    kap  = dat.sampled_segments{1,i}.kappa;
    cr=[cr;kap];
    
    cc = dat.sampled_segments{1,i}.chain;
    chain_number = [chain_number,cc];
    
    ll = (dat.sampled_segments{1,i}.ind) -1;    
    l_on_chain = [l_on_chain;ll];
end

curv_matrix = [cr,chain_number',l_on_chain];
savetxt_v01(filename,curv_matrix,'kappa');

f = figure();
hist(cr,100)
set(gca,'fontsize',14)
xlabel('Kappa (nm^{-1})')
ylabel('#')
saveimage(f,'kappa','eps',FN);

f = figure();
hist(1./cr,100)
set(gca,'fontsize',14)
xlabel('Radius of Curvature (nm)')
ylabel('#')
saveimage(f,'curvature','eps',FN);

close all

