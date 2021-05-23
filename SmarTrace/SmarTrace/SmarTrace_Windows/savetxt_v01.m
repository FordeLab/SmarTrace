function savetxt_v01(filename,data,dataname)

FN = strrep(filename(1:end-4),'_','-');
FN_path = [pwd '/' FN '/' '%s.txt'];

dlmwrite(sprintf(FN_path,dataname),data);
