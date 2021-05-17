
function SmTr_ChainSamplingMOD(FN,ex)
% SmTr_analysis1 This code samples the traced chains from SmarTrace tracking 
% program using a modiefied version of Faas et al (2009)
% sampling method. It randomly samples chain pieces with lengths
% drawn from a predefined set such that the whole chain is used without 
% overlapping. The random sizes of chain segments are shuffled to avoid 
% concentration of smaller lengths at other end of the chain. 
% Data saved from this step can be used with with SmTr_Analysis.m for 
% statistical analysis of the chains.
% copyright 2016 Naghmeh Rezaei
% *******Do NOT distribute. *******
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Changes by Mathew Schneider, 2018
% added if statement to for loop to check for empty structures
% should now work if chains have been missed
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
method = 'randperm'; %choose among:  'randperm', 'dropfirstlast', 'randstartend'
rand_draws = 50;
bin_size = 10;
max_size = 200;
outputstruct = struct();
outputstruct.sampled_segments = cell(0,0);

pixel_trim = 5; %first/last pixels to be excluded from the chains

data = load(FN);
n_chains = 0;
if ex ~= 0
    data.bsplines_norm(ex) = [];
    data.der(ex) = [];
    
end

bspline_norm = struct2cell(data.bsplines_norm);
bspline_norm = bspline_norm(~cellfun('isempty',bspline_norm));
tempDer = struct2cell(data.der);
for i=1:size(tempDer,3)
    Xder{i} = tempDer{1,1,i};
    Yder{i} = tempDer{2,1,i};
end
Xder = Xder(~cellfun('isempty',Xder));
Yder = Yder(~cellfun('isempty',Yder));
conLengths = struct2cell(data.contour_lengths);
conLengths = cell2mat(conLengths);
num_chains = length(bspline_norm);

RandIndex = randperm(num_chains);

for i = 1:(num_chains)
    display(['Chain no ' num2str(i) ' / ' num2str(num_chains)]);   
    
        ff = RandIndex(i);
    
        x_norm = bspline_norm{1,1,ff}(1,pixel_trim:end-pixel_trim)';
        y_norm = bspline_norm{1,1,ff}(2,pixel_trim:end-pixel_trim)';
        
    
    
        x_der = Xder{1,ff}(pixel_trim:end-pixel_trim);
        y_der = Yder{1,ff}(pixel_trim:end-pixel_trim);

        dir_der = [x_der,y_der];


        dir_der2 = dir_der;
        dir_der2(:,3) = 0;
        dir_derc = cross(dir_der2(2:end,:), dir_der2(1:end-1,:));
        dir_derd = dot(dir_der2(2:end,:)', dir_der2(1:end-1,:)')';
        % cumulative angle differences
        dir_cumder = cumsum(atan2(dir_derc(:,3), dir_derd));
        dir_cumder = [0 ; dir_cumder];

        %arclen = arclength(x_norm, y_norm);
        arclen= length(x_norm)-1;

        R2s_mean = [];
        R2s_std = [];

        cosines = [];
        seps = [];
        R2s = [];

        thetas = [];
        seps_angles = [];

        % separations = 1:1:max(arclen);
        for nn = 1:rand_draws

            sz=[];

            while length(sz)==0
                [ind, sz] = sampling(arclen, [bin_size:bin_size:max_size], method);
            end


            if length(ind)-length(sz) ~= 1
                error('here');
            end        

            x_points = x_norm(ind);
            y_points = y_norm(ind);


            theta = diff(dir_cumder(ind)); 

            thetas = [thetas;theta];
            seps_angles = [seps_angles;sz'];

            %cosine = dot(der(1:end-1,:),der(2:end,:),2);
            cosine = cos(theta);
            cosines = [cosines;cosine];
            seps = [seps;sz'];


            %end-to-end (R^2) calculations
            rr = [x_points,y_points];
            R2 = sum(diff(rr).^2,2);
            R2s = [R2s;R2];

            for ii=1:length(sz)
                n=length(outputstruct.sampled_segments)+1;
                outputstruct.sampled_segments{n}.chain = ff;
                outputstruct.sampled_segments{n}.draw = nn;
                outputstruct.sampled_segments{n}.sep = sz(ii);
                outputstruct.sampled_segments{n}.ind = ind(ii);
                outputstruct.sampled_segments{n}.theta = theta(ii);
                outputstruct.sampled_segments{n}.cosine = cosine(ii);
                outputstruct.sampled_segments{n}.R2 = R2(ii);
            end
        end


        corel_single = [cosines, seps];
        thetas_single = [thetas,seps_angles];

        end2cont_single = [R2s,seps];
        goodcontour = arclen;

        total_end2end = sqrt (( y_norm(end) - y_norm(1) ).^2 +...
            (( x_norm(end) - x_norm(1) ).^2 ));

        add_to_structure(ff,corel_single,end2cont_single,goodcontour,...
            total_end2end,thetas_single,n_chains);
    

end

save_data(FN);


    function add_to_structure(n,corel_single,end2cont_single,goodcontour...
            ,total_end2end,thetas_single,n_chains)
        
        outputstruct.contour_lengths(n).contourL = goodcontour;
        
        outputstruct.mat_correlations(1,n).corelfib = corel_single;
        outputstruct.mat_wormlike(n).wormfib = end2cont_single;
        
        outputstruct.length_conformation(n, 1) = goodcontour;
        outputstruct.length_conformation(n, 2) = total_end2end;
        
        outputstruct.n_chains = n_chains;
        
        %outputstruct.thetas(n).thetas = thetas_single;
    end


    function save_data(FN)
        save(['Re1-' method '-' num2str(rand_draws) '_' num2str(bin_size) ...
            '_' num2str(max_size) '_' FN], '-struct', 'outputstruct');
    end
end


