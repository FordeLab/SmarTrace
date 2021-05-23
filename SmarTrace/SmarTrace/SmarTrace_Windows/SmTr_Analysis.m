
function SmTr_Analysis(tempfilename)
% SmTr_analysis2 This code provides statistical analysis of the sampled 
% chains from SmTr_ChainSampling.m.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Edits by Mathew Schneider:
%Line 314 - changed the weighting to 1/sem^2
%Line 233 - changed weighting to 1/w^2
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% clear
% clc

% in current folder
file_name = tempfilename;
% number of bins for angle in each dl
nbin_angles = 50;
% bin size for lengths
dl_lengths = 10;
% WLC fit lower limit
low_limit = 10;
% WLC fit upper limit
high_limit = 200;

%angle range over which the histogram is calculated
rng = pi; %choose either pi/2 or pi

% TODO: make these functions
% % bin number for each dl segments
% length_to_bin = @(l) floor(l/dl_lengths)+1;
% % returns the centre of a given bin
% bin_to_length = @(b) (b-0.5)*dl_lengths;
% bin number for each dl segments
length_to_bin = @(l) ceil(l/dl_lengths);
% returns the centre of a given bin
bin_to_length = @(b) b*dl_lengths;
% fixed bin sizes
angle_to_bin = @(a) floor(a*nbin_angles/(pi+0.0001))+1;
% returns the angle given the bin number
bin_to_angle = @(b) (b-0.5)*(pi+0.0001)/nbin_angles;


% angle bining both with specified range
nbin_cumangles = 50;
cumangle_to_bin = @(a) floor(a*nbin_angles/(pi+0.0001))+1;


cos_nbins = 20;
cos_to_bin = @(x) ceil((x+1)*cos_nbins/2);
bin_to_cos = @(b) (b*2/cos_nbins)-1;

FN = strrep(file_name(1:end-4),'_','-');
% FN = (strcat('v7-',strrep(file_name(1:end-8),'_','-chains')))%,'-dl'),sprintf('%d',dl_lengths)));

dirname = FN;
mkdir(dirname);

% curv_v01(file_name);
run();

    function run()
        
        sampled_segments = read_reanalyze_data(file_name);
        
%         angle_sep_2d_hist = accumcells(sampled_segments, ...
%             @(d) [length_to_bin(d.sep); angle_to_bin(d.theta)], ...
%             'theta', ...
%             @(x) length(x), ...
%             0);
        

        sep_bin_count = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'theta', @(x) length(cos(x)), NaN);
        
        cor_binned_mean = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'theta', @(x) mean(cos(x)), NaN);
        cor_binned_std = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'theta', @(x) std(cos(x)/sqrt(length(x)),1), NaN);
%         /sqrt(length(x))
        r2_binned_mean = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'R2', @mean, NaN);
        r2_binned_std = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'R2', @(x) std(x,1)/sqrt(length(x)), NaN);
        
        theta2_binned_mean = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'theta', @(x) mean(x.^2), NaN);
        theta2_binned_std = accumcells(sampled_segments, @(d) ...
            length_to_bin(d.sep), 'theta', @(x) std(x.^2,1)/sqrt(length(x)), NaN);
        
%         cor_bin_2D = accumcells(sampled_segments, @(d) ...
%             [length_to_bin(d.sep); d.cosine], 'theta', @(x) length(x), 0);
%         R2_bin_2D = accumcells(sampled_segments, @(d) ...
%             [length_to_bin(d.sep); d.R2], 'theta', @(x) length(x), 0);
        
        plot_errordata('<cos>', cor_binned_mean, cor_binned_std);
        plot_errordata('<R2>', r2_binned_mean, r2_binned_std);
        plot_errordata('<theta2>', theta2_binned_mean, theta2_binned_std);
        
        
        %WLC fits
        [lp,ef,gof]=fit_WLC('<cos>', bin_to_length(1:length(cor_binned_mean))',...
            cor_binned_mean, cor_binned_std, sep_bin_count,low_limit, high_limit,[],1)
        [lp,ef,gof]=fit_WLC('<R2>', bin_to_length(1:length(r2_binned_mean))'...
            , r2_binned_mean, r2_binned_std, sep_bin_count,low_limit, high_limit,[],1)
        fit_WLC('<theta2>', bin_to_length(1:length(theta2_binned_mean))'...
            , theta2_binned_mean, theta2_binned_std,sep_bin_count,low_limit, high_limit,[],1);
%         
        
        WLCfit_results = [];
        
        dr = 1;
        for fit_range = 5*dl_lengths:dl_lengths:high_limit
            row = [fit_range];
            [lp,ef,gof] = fit_WLC('<cos>', bin_to_length(1:length(cor_binned_mean))',...
                cor_binned_mean, cor_binned_std,sep_bin_count,low_limit, fit_range,[], dr)
            row = [row lp ef gof];
            [lp,ef,gof] = fit_WLC('<R2>', bin_to_length(1:length(r2_binned_mean))'...
                , r2_binned_mean, r2_binned_std,sep_bin_count,low_limit, fit_range,[], dr)
            row = [row lp ef gof];
            [lp,ef,gof] = fit_WLC('<theta2>', bin_to_length(1:length(theta2_binned_mean))'...
                , theta2_binned_mean, theta2_binned_std, sep_bin_count,low_limit, fit_range,[], dr);
            row = [row lp ef gof];
            WLCfit_results = [WLCfit_results; row];
        end
        close all
        savetxt({'range' 'coslp' 'cosef' 'cosgof' 'r2lp' 'r2ef' 'r2gof' 't2lp' 't2ef' 't2gof'}, WLCfit_results, 'fit_range');
        
        f = figure();
        errorbar(WLCfit_results(:,1),WLCfit_results(:,2),WLCfit_results(:,3),'bo');
        hold on
        errorbar(WLCfit_results(:,1),WLCfit_results(:,5),WLCfit_results(:,6),'ro');
%         errorbar(WLCfit_results(:,1),WLCfit_results(:,8),WLCfit_results(:,9),'ko');
        set(gca,'fontsize',18);
        xlabel( 'Fit Range (nm)' );
        ylabel( 'Persistence Length (nm)' );
        legend('<cos>', '<R^2>')
        t = 'fit_range'; 
        title('Recombinant Yeast Collagen I')
        saveimage(f,t,'eps',dirname);
        hold off
%     
calculate_angle_hist(sampled_segments, rng );
        
        % histogram of sample segment lengths
        seg_hist = accumcells(sampled_segments, @(d) [length_to_bin(d.sep)], 'sep', @(x) length(x), 0);
        X = bin_to_length(1:size(seg_hist,1));
        
        f= figure();
        h = loglog(X,seg_hist,'ko');
        set(h,'markersize',8)
        set(gca,'fontsize',18);
        xlabel( 'Segment length (nm)' );
        ylabel( '#' );
% %         title(['seg_histt-' FN]);
        t = 'seg_hist';
        saveimage(f,t,'eps',dirname);
       
%         calculate_fit_Lmax(cor_binned_mean, r2_binned_mean, theta2_binned_mean);
%                 kurt = accumcells(sampled_segments, @(d) [d.draw; length_to_bin(d.sep)], 'theta', @(x) mean(x.^4)/(mean(x.^2))^2, NaN);
%         plot_kurtosis(kurt);
                
%         seg_along_chain = accumcells(sampled_segments, @(d) [length_to_bin(d.sep)], 'ind', @mean, NaN);
%         X = bin_to_length(1:size(seg_along_chain,1))
%         Y = seg_along_chain' + X/2;
%         plot(X,Y);
        
        
    end

    function sampled_segments = read_reanalyze_data(file_name)
        % read all segment data from reanalyze code
        data_struc = load(file_name);
        sampled_segments = data_struc.sampled_segments;
    end


    function angle_hist = calculate_angle_hist(segs, rng )        
        ind = [2,3,4,5,6,7]; 
% ind = 1;
        len_bin = ind*dl_lengths;
        angle_hist = accumcells(segs, @(d) [length_to_bin(d.sep); ...
            angle_to_bin(abs(mod(d.theta+rng,2*rng)-rng))], 'theta', @(x) length(x), 0);
        angle_hist = angle_hist';
        sep_bin_count = sqrt(angle_hist);
        %sum of all angle probability in each bin length
        sum_prob = repmat(sum(angle_hist,1),size(angle_hist,1),1);
        
        x_range = (bin_to_angle(1:nbin_angles))';
        
        angle_prob = angle_hist./sum_prob;
        lnG = -log(angle_prob);
        
        plot(x_range(1:size(angle_hist,1)), lnG(:,ind), 'o');               
        row = [];
        lnGfit_results = [];
        for i = 1:length(ind)
            row = [len_bin(i)];
            y = lnG(:,ind(i));
            y = y(~isinf(y));
            w = sep_bin_count(:,ind(i));
            w = w(~isinf(y));
            
            [lp,ef,gof] = fit_lnG('lnG',x_range(~isinf(y)), y,w, len_bin(i), 1);
            row = [row lp ef gof];
            lnGfit_results = [lnGfit_results; row];
        end
        
        close all
        
        plot(1./lnGfit_results(:,1),lnGfit_results(:,3));
        fit_lnG('lnG_lp',1./lnGfit_results(:,1),lnGfit_results(:,3),[],[],1)
        
        regexp(sprintf(['L%f '],[2,4,5]*2.15), ' ', 'split')
        
        FN_path1 = [pwd '/' dirname '/' 'angle_prob_%s.txt'];
        
%        dlmwrite(sprintf(FN_path1,'angle_prob'), [x_range(1:size(angle_hist,1)) ...
%            angle_prob]);
        
%        dlmwrite(sprintf(FN_path1,'lnG'), [x_range(1:size(angle_hist,1)), ...
%            lnG]);
    end

    function [B,ef,gof] = fit_lnG(datatype,x, y, w, l, plt)
       
        if strcmp(datatype,'lnG')
            
            ft = fittype( 'A+(B*x^2)', 'independent', 'x', 'dependent', 'y');
            tit = datatype;
            
            opts = fitoptions(ft);
            opts.Display = 'Off';
            %         opts.StartPoint = 20;
            w = 1./w.^2;
            w(isinf(w))=NaN;
            
            y_label = sprintf('-ln(P (\\theta, %d nm))', l);
            x_label = 'Angle (Rad)';
            
        elseif strcmp(datatype,'lnG_lp')
            
            ft = fittype( '0.5*p*x', 'independent', 'x', 'dependent', 'y');
            tit = datatype;
            
            y_label = sprintf('-ln(P) fit parameter');
            x_label = 'Length^{-1} (nm^{-1})';
            
            opts = fitoptions(ft);
            opts.Display = 'Off';
            w = [];

        end
            
        [fitresults, gof, output] = fit(x, y, ft, opts);
        B = coeffvalues(fitresults);
        gof = gof.rsquare;
        ci = confint(fitresults);
        ef = (ci(2)-ci(1))/2;
    
        if plt
            f = figure();
            h = plot( fitresults, x, y,'ko');
            hold on;
            if strcmp(datatype,'lnG')
                errorbar(x, y, w, 'ko');
            end

            set(gca,'fontsize',18);
            xlabel( x_label );
            ylabel( y_label );
            %         annotation('textbox',[0.2, 0.8, .1,.1], 'String', corr_fit);
            saveimage(f, ['WLC-' tit '_' num2str(l)], '-depsc',dirname);
            
            figure; plot(output.residuals,'o')
            set(gca,'fontsize',18,'fontweight','bold');
            xlabel( 'Length (nm)' );
            ylabel( y_label );
        end
    end

    %function [lp,rw,Jw,Sigmaw,msew] = fit_WLC(datatype, x, y, w, low, high, l, plt)
    function [lp,ef,gof] = fit_WLC(datatype, x, y, sem,w, low, high, l, plt)
        

        % fits correlation and R^2 data with no bootstrapping
        
        if strcmp(datatype, '<cos>')
            ft = fittype( 'exp(-x/(2*p))', 'independent', 'x', 'dependent', 'y');
            y_label = '<cos \theta>';
            tit = 'cos';
        elseif strcmp(datatype, '<R2>')
            ft = fittype( '4*p*x*(1-(2*p/x)*(1-exp(-x/(2*p))))', 'independent',...
                'x', 'dependent', 'y');
            y_label = '<R^2> (nm^2)';
            tit = 'R2';
        else
            ft = fittype('x/p', 'independent',...
                'x', 'dependent', 'y');
            y_label = '<\theta^2>';
            tit = 'theta2';
        end
        
        opts = fitoptions(ft);
        opts.Display = 'Off';
        opts.Lower = low;
        opts.StartPoint = 20;
        opts.Upper = high;
        ex = excludedata( x, y, 'domain', [low  high] );
        opts.Exclude = ex;
% x = x(~ex);
% y = y(~ex);
% w = w(~ex);
% sem = sem(~ex);

%         opts.Weights = w / mean(w);
opts.Weights = 1./sem.^2;
       
        [fitresults, gof, output] = fit( x, y, ft, opts );
        lp = coeffvalues(fitresults);
        gof = gof.rsquare;
        ci = confint(fitresults);
        ef = (ci(2)-ci(1))/2;
        
        if plt
            f = figure();
            if max(ex)
                h = plot( fitresults, x, y,ex);
               
            else
                h = plot( fitresults, x, y);
                set(h, 'MarkerSize', 8);
                set(h, 'linewidth', 2);
                set(h(1),'Marker','o');

                set(h(1),'Markeredgecolor','k');
            end
            hold on;
            
            if max(ex)
                errorbar(x(ex), y(ex), sem(ex), 'g+');
            end
            h2= errorbar(x(~ex), y(~ex), sem(~ex), 'k');
            set(h2(1),'Marker','+');
            set(h(1),'marker','o');
            set(h2(1),'Linestyle','none');

            set(gca,'fontsize',20,'fontweight','bold');
            xlabel( 'Length (nm)' );
            ylabel( y_label );
            xlim([0 high_limit]) 
%             if strcmp(datatype,'<R2>')
%                 leg = findobj (gcf, 'tag', 'legend');
%                 set(leg,'location','northwest')                
%             end
%             title(['WLC fit-' FN])
% legend('Experiment','WLC fit')
            hold off
            %         annotation('textbox',[0.2, 0.8, .1,.1], 'String', corr_fit);
            saveimage(f, ['WLC-' tit '-' num2str(high)], '-depsc',dirname);
            
            close all
            f= figure; errorbar(x(~ex),output.residuals,sem(~ex),'o')

            set(gca,'fontsize',18);
            xlabel( 'Length*10 (nm)' );
            ylabel( 'Residuals, yFit - y' );
            saveimage(f, ['Residual-' tit '-' num2str(high)], '-depsc',dirname);

        end
        
        return
        
        % % new method
        
        % fits correlation and R^2 data
        
%         w = w / mean(w);
%         yw = sqrt(w).*y;
        
        start = 20;

        if strcmp(datatype, '<cos>')
            
            modelFun = @(p,x) exp(-x/(2*p));
            modelFunw = @(p,x) sqrt(w).*modelFun(p,x);            
            
            y_label = '<cos \theta>';
            x_label = 'Length (nm)';
            tit = 'cos';
            
        elseif strcmp(datatype, '<R2>')
            
            modelFun = @(p,x) 4*p*x*(1-2*p/x*(1-exp(-x/2/p)));
            modelFunw = @(p,x) sqrt(w).*modelFun(p,x);            
            
            y_label = '<R^2>';
            x_label = 'Length (nm)';
            tit = 'R2';

        elseif strcmp(datatype, '<theta2>')
            
            modelFun = @(p,x) x/p;
            modelFunw = @(p,x) sqrt(w).*modelFun(p,x);
            
            y_label = '<\theta^2>';
            x_label = 'Length (nm)';
            tit = 'theta2';           
            
        end
        
        [lp,rw,Jw,Sigmaw,msew] = nlinfit(x,yw,modelFunw,start);
        lp
        
        rmsew = sqrt(msew);
        bCIw = nlparci(lp,rw,'cov',Sigmaw);
        seFitw = sqrt(diag(Sigmaw));
        
        xgrid = linspace(min(x),max(x),100)';
        [yFitw, deltaw] = nlpredci(modelFun,xgrid,lp,rw,'cov',Sigmaw); 
       
        if plt
            
            f = figure();
            h = plot(x,y,'bo', xgrid,yFitw,'r-',xgrid,yFitw+deltaw,...
                'b:',xgrid,yFitw-deltaw,'b:');
            legend({'Data', 'WLC fit', '95% Confidence Limits'},...
                'location','NorthEast');
            
            set(gca,'fontsize',18);            
            set(h, 'MarkerSize', 7);
            set(h(1),'Marker','o');
            
            xlabel(x_label); ylabel(y_label);
            
            saveimage(f,['WLC-' tit],'-depsc',dirname);
            
            %compare weighted vs unweighted fit
            f = figure();
            
            [bFit,r,J,Sigma,mse] = nlinfit(x,y,modelFun,start);
            [yFit,delta] = nlpredci(modelFun,xgrid,bFit,r,'cov',Sigma);
            
            h= plot(x,y,'ko', xgrid,yFitw,'b-',xgrid,yFit,'r-');
            legend({'Data', 'Weighted fit', 'Unweighted fit'},'location','NorthEast');
            set(gca,'fontsize',18);
            set(h, 'MarkerSize', 7);
            xlabel(x_label); ylabel(y_label);
            
            saveimage(f,['cmpr_weight_fit' tit], '-depsc',dirname);
            
            %residual analysis
            f = figure();
            h = plot(x,rw.*sqrt(w),'b^');
            if feature('HGUsingMATLABClasses')
                hy = specgraphhelper('createConstantLineUsingMATLABClasses',...
                    'LineStyle',':','Color',[.5 .5 .5],'Parent',gca);
                hy.Value = 0;
            else
                graph2d.constantline(0,'linestyle',':','color',[.5 .5 .5]);
            end
            set(gca,'fontsize',18);
            set(h, 'MarkerSize', 7);
            xlabel('x'); ylabel('Residuals, yFit - y');
            set(gca,'fontsize',18);
            set(h, 'MarkerSize', 7);
            saveimage(f,['Residual' tit],'-depsc',dirname);
            
        end
        
    end


    function plot_errordata(datatype, ymean, ystd)

        if strcmp(datatype, '<cos>')
            y_label = '<cos(\theta)>';
            tit = 'cos';
        elseif strcmp(datatype, '<R2>')
            y_label = '<R^2> (nm^2)';
            tit = 'R2';
        else
            y_label = '<\theta^2> (Rad^2)';
            tit = 'theta2';
        end
        
        f = figure();
        data = [bin_to_length(1:length(ymean))', ymean, ystd];
        errorbar(data(:,1), data(:,2), data(:,3), 'o');
        set(gca,'fontsize',18);
        xlabel('Length (nm)');
        ylabel(y_label)
        
%         title(FN);
        saveimage(f,tit,'esp',dirname);
        
        %save the avg binned results
        FN_path1 = [pwd '/' dirname '/' '_Results_%s.txt'];
%        dlmwrite(sprintf(FN_path1,datatype),data);
    end

    function plot_angle_hist(angle_hist)
        
        f = figure();
        hold on
        
        for ldraw = [low_limit 2*low_limit 3*low_limit 4*low_limit];
            m = nanmean(angle_hist, 1);
            s = nanstd(angle_hist, 1);
            
            %normalize = repmat(sum(hist_count,2),1,size(hist_count,2));
            
            errorbar((bin_to_angle(1:nbin_angles))*180/pi, m(1,length_to_bin(ldraw),:), s(1,length_to_bin(ldraw),:));
        end
        xlabel('Angle (Degrees)');
        ylabel('Count');
        
        saveimage(f,'angle_dist','eps',dirname);
        
        end

    function plot_kurtosis(kurt)
        %plot kurtosis of angles as a function of length
        
        f = figure();
        plot(bin_to_length(1:size(kurt,2)), nanmean(kurt, 1),'o');
        set(gca,'fontsize',18);
        xlabel('Length (nm)');
        ylabel('Kurtosis');
        
        hold on;
        plot(bin_to_length([1 size(kurt,2)]),[3 3], 'r-')
        plot(bin_to_length([1 size(kurt,2)]),[1.8 1.8], 'g-')
        
        %         axis([0 160 0 10]);
        title(FN);
        hold off
        saveimage(f,'kurtosis','eps',dirname);
    end

    function calculate_fit_Lmax(cor, r2,theta2)
        
        lp_cor = [];
        gof_cor = [];
        lp_r2 = [];
        gof_r2 = [];
        L_max = [];
        lp_theta2 = [];
        gof_theta2 = [];
        
        cor_len = bin_to_length(1:size(cor,2));
        r2_len = bin_to_length(1:size(r2,2));
        theta2_len = bin_to_length(1:size(theta2,2));
        
        ll = 0;
        for lmax = 50:1*dl_lengths:high_limit
            ll = ll+1;
            L_max(ll,1) = lmax;
            
            for draw = 1:size(cor,1)
                
                [lp_cor(draw,ll), gof_cor(draw,ll)] = fit_WLC('<cos>',cor_len, cor(draw,:), low_limit, lmax,0);
                [lp_r2(draw,ll), gof_r2(draw,ll)] = fit_WLC('<R2>', r2_len, r2(draw,:), low_limit, lmax,0);
                %                 [lp_theta2(draw,ll), gof_theta2(draw,ll)] = fit_WLC('<theta2>', theta2_len, theta2(draw,:), low_limit, lmax,0);
                
            end
        end
        
        f = figure();
        errorbar(L_max,nanmean(lp_cor,1),nanstd(lp_cor,1,1),'bo');
        hold on;
        errorbar(L_max,nanmean(lp_r2,1),nanstd(lp_r2,1,1),'ro');
        %         errorbar(L_max,nanmean(lp_theta2,1),nanstd(lp_theta2,1,1),'go');
        legend('<cos(\theta)>' , '<R^2>')%, '<\theta^2>');
        xlabel('Lmax (nm)');
        ylabel('Persistence Length (nm)');
        grid on
        hold off
        saveimage(f,'fits_Lmax','eps',dirname);
        
    end
 
    function savetxt(headers, data, name)
        FN_path = [pwd '\' dirname '\' FN '_' name '.txt'];
        fileID = fopen(FN_path,'w');
        for s=headers
            fprintf(fileID, '%s ',s{1});
        end
        fprintf(fileID, '\n');
        fprintf(fileID, [repmat('%f ',1, length(headers)), '\n'], data');
        fclose(fileID);
    end

end



