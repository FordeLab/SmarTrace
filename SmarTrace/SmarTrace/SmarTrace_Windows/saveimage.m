
function saveimage(fig,name,format,FN)

dirname = FN;
% mkdir(dirname);
        
        if strcmp(format,'fig')        
        FN_path = [pwd '\' dirname '\' name '.fig'];        
        saveas(fig,FN_path,'fig');
        
        else
            FN_path = [pwd '\' dirname '\' name '.png'];
%             print(fig,FN_path,'-depsc');
 print(fig,FN_path,'-dpng');
        end