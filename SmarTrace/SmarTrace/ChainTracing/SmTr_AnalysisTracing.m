
function P = SmTr_AnalysisTracing(SampledStruct,LL,ko)
% SmTr_analysis2 This code provides statistical analysis of the sampled 
% chains from SmTr_ChainSampling.m.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Edits by Mathew Schneider (2018):
%Line 314 - changed the weighting to 1/sem^2
%Line 233 - changed weighting to 1/w^2
%General - added ability to do curved chain fits using the cWLC model
%          set the ko input to 0 for WLC fits, 1 for cWLC fits
%Saving - fewer images are saved to reduce runtime and space
%         only WLC fits, data, and residuals for 200 nm are saved
%         also changed the data format to .csv from .txt
%GOF - now uses reduced chi squared for gof instead of r squared and also 
%      calculates the Akaike Information Criterion (AIC)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% clear
% clc



% number of bins for angle in each dl
nbin_angles = 50;
% bin size for lengths
dl_lengths = 10;
% WLC fit lower limit
low_limit = LL;
% WLC fit upper limit
high_limit = 200;
% WLC fit start point
start_point = 80;
% cWLC lower limit for k then lp
curve_low = [0 0];
% cWLC high limit for k then lp
curve_high = [.03 1000];
% cWLC start points for k then lp
curve_start = [0 80];

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



% curv_v01(file_name);
% clc
       

        sampled_segments = SampledStruct.sampled_segments;
        
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
        
        
        %WLC fits
        [P,ef,gof]=fit_WLC('<cos>', bin_to_length(1:length(cor_binned_mean))',...
            cor_binned_mean, cor_binned_std, sep_bin_count,low_limit, high_limit,[],1);
        [P,ef,gof]=fit_WLC('<R2>', bin_to_length(1:length(r2_binned_mean))'...
            , r2_binned_mean, r2_binned_std, sep_bin_count,low_limit, high_limit,[],1);



        %function [lp,rw,Jw,Sigmaw,msew] = fit_WLC(datatype, x, y, w, low, high, l, plt)
    function [lp,ef,gof] = fit_WLC(datatype, x, y, sem,w, low, high, l, plt)
        

        % fits correlation and R^2 data with no bootstrapping
    if ko == 0    
        
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
        opts.Lower = 0;
        opts.StartPoint = start_point;
        opts.Upper = 1000;
        ex = excludedata( x, y, 'domain', [low  high] );
        opts.Exclude = ex;
    elseif ko == 1
        if strcmp(datatype, '<cos>')
            ft = fittype( 'cos(ko*x)*exp(-x/(2*p))',...
                'independent', 'x', 'dependent', 'y');
            y_label = '<cos \theta>';
            tit = 'cos';
        elseif strcmp(datatype, '<R2>')
            ft = fittype( ['(4*x*p/(1+4*ko^2*p^2)^2)*',...
                '(1-2*(p/x)*(1-4*ko^2*p^2)',...
                '*(1-cos(ko*x)*exp(-x/2/p))+(4*ko*p^2/x)',...
                '*(ko*x -2*sin(ko*x)*exp(-x/2/p)))'], 'independent',...
                'x', 'dependent', 'y');
            y_label = '<R^2> (nm^2)';
            tit = 'R2';
        else
            ft = fittype('x*(ko^2*x*p+1)/p', 'independent',...
                'x', 'dependent', 'y');
            y_label = '<\theta^2>';
            tit = 'theta2';
        end
        
        opts = fitoptions(ft);
        opts.Display = 'Off';
        opts.Lower = curve_low;
        opts.StartPoint = curve_start;
        opts.Upper = curve_high;
        ex = excludedata( x, y, 'domain', [low  high]);
        opts.Exclude = ex;
    end
    
% x = x(~ex);
% y = y(~ex);
% w = w(~ex);
% sem = sem(~ex);

%         opts.Weights = w / mean(w);
opts.Weights = 1./sem.^2;
       
        [fitresults, gof, output] = fit( x, y, ft, opts );
        lp = coeffvalues(fitresults);
        
        if ko == 1
        gof = gof.sse./18;
        elseif ko == 0
            gof = gof.sse./19;
        end
        
        ci = confint(fitresults);
        ef = (ci(2,:)-ci(1,:))/2;
        
        
        
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

        
    end





end
