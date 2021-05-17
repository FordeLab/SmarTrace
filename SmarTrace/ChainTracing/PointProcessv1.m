function chains = PointProcess(input_image)
% 
% %Load in image to label
% Im = rgb2gray(imread(input_image));%Test image
% Im = Im(19:530,38:549);
Im = input_image;

%Load network and label image
netv5 = load('NetworkGen2.mat');%Better labeling and changed learning rate
I = semanticseg(Im,netv5.net);% Classify pixels with network
B = labeloverlay(Im,I);% overlay classification with image

%Turn the categorical label matrix into a binary image
%   then remove small and large objects
bw = logical(ismember(I,'Chains')+ismember(I,'Overlap'));
 for n = 1:3
    for i=2:length(bw)-1
        for j=2:length(bw)-1
            ul = bw(i,j-1) + bw(i+1,j) + bw(i+1,j-1);
            ur = bw(i,j+1) + bw(i+1,j) + bw(i+1,j+1);
            bl = bw(i,j-1) + bw(i-1,j) + bw(i-1,j-1);
            br = bw(i,j+1) + bw(i-1,j) + bw(i-1,j+1);
            if ul < 3 && ur < 3 && bl < 3 &&  br < 3
                bw(i,j) = 0;
            end
        end
    end
 end
bw = bwareaopen(bw,80);
chains = regionprops(bw,'Area');
cc = bwconncomp(bw,8);
L = labelmatrix(cc);
Lbw = ismember(L,find([chains.Area] > 500));
bw = ismember(L,find([chains.Area] <= 500));
chains = regionprops(bw,'Area','BoundingBox');
Lbw = bwmorph(Lbw,'open');

%Remove pixels with less than 3 adjacent white pixels
 for n = 1:3
    for i=2:length(bw)-1
        for j=2:length(bw)-1
            ul = bw(i,j-1) + bw(i+1,j) + bw(i+1,j-1);
            ur = bw(i,j+1) + bw(i+1,j) + bw(i+1,j+1);
            bl = bw(i,j-1) + bw(i-1,j) + bw(i-1,j-1);
            br = bw(i,j+1) + bw(i-1,j) + bw(i-1,j+1);
            if ul < 3 && ur < 3 && bl < 3 &&  br < 3
                bw(i,j) = 0;
            end
        end
    end
 end

%Remove small objects again and get object properties
Circ = regionprops(bw,'Circularity','FilledImage','BoundingBox');

for i = 1:length(Circ)
    box = floor(Circ(i).BoundingBox);
    if Circ(i).Circularity > 0.4
        bw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = ...
            bw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) - Circ(i).FilledImage;
    end
end
bw = bwareaopen(bw,150,4);
chains = regionprops(bw,'BoundingBox','FilledImage');

for i = 1:length(chains)
    skel = bwskel(chains(i).FilledImage,'MinBranchLength',30);
    ep = bwmorph(skel,'branchpoints');
    [r c] = find(ep);
    if ~isempty(r)
        for j = 1:length(r)
            temp = zeros(size(chains(i).FilledImage)+13);
            temp(7:end-7,7:end-7) = chains(i).FilledImage;
            temp(r:r+10,c:c+10) = 0;
            chains(i).FilledImage = temp(7:end-7,7:end-7);
        end
        box = floor(chains(i).BoundingBox);
        bw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = chains(i).FilledImage;
    end
end

Lchains = regionprops(Lbw,'Area','BoundingBox','FilledImage');
for i = 1:length(Lchains)
    skel = bwskel(Lchains(i).FilledImage,'MinBranchLength',10);
    ep = bwmorph(skel,'branchpoints');
    [r c] = find(ep);
    if ~isempty(r)
        for j = 1:length(r)
            temp = zeros(size(Lchains(i).FilledImage)+13);
            temp(7:end-7,7:end-7) = Lchains(i).FilledImage;
            temp(r:r+10,c:c+10) = 0;
            Lchains(i).FilledImage = temp(7:end-7,7:end-7);
        end
        box = floor(Lchains(i).BoundingBox);
        Lbw((box(2)+1:box(2)+box(4)),(box(1)+1:box(1)+box(3))) = Lchains(i).FilledImage;
    end
end
bw = bwareaopen(bw,100,8);
Lbw = bwareaopen(Lbw,100,8);
bw = logical(bw+Lbw);

chains = regionprops(bw,'BoundingBox','FilledImage');

for i = 1:length(chains)
    [r1 c1] = size(chains(i).FilledImage);
    pic = chains(i).FilledImage;
    box = chains(i).BoundingBox;

    %Get the skeleton of the chain
    skel = bwskel(pic,'MinBranchLength',30);    
    chains(i).skeleton = skel;
    [r c] = find(skel);
    N = length(r);
    bp = bwmorph(skel,'branchpoints');
    if all(bp(:)==0) && length(r) > 30
        for j = 1:N/2
                [r c] = find(skel);
                endpnts = bwmorph(skel,'endpoints');
                if j == 1
                    [y x] = find(endpnts);    
                    chains(i).points(1,j) = x(1);
                    chains(i).points(2,j) = y(1);
                    chains(i).points(1,N-j) = x(2);
                    chains(i).points(2,N-j) = y(2);
                    skel = skel - endpnts;
                elseif j > 1
                    [y x] = find(endpnts);

                    if abs(x(1) - chains(i).points(1,j-1)) > 2 || ...
                            abs(y(1) - chains(i).points(2,j-1)) > 2 ...
                            abs(x(2) - chains(i).points(1,N-j+1)) > 2 ...
                            || abs(y(2) - chains(i).points(2,N-j+1)) > 2;

                        chains(i).points(1,j) = x(2);
                        chains(i).points(2,j) = y(2);
                        chains(i).points(1,N-j) = x(1);
                        chains(i).points(2,N-j) = y(1);
                        skel = skel - endpnts;
                    else
                        chains(i).points(1,j) = x(1);
                        chains(i).points(2,j) = y(1);
                        chains(i).points(1,N-j) = x(2);
                        chains(i).points(2,N-j) = y(2);
                        skel = skel - endpnts;
                    end
                end
        end
        box = floor(box);
        chains(i).points(1,:) = chains(i).points(1,:) + box(1);
        chains(i).points(2,:) = chains(i).points(2,:) + box(2);
    else
        chains(i).points = [];
    end
end

end

