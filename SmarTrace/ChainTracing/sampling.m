% cut track into pieces of random "length"
% METHOD is either 'pla', 'randperm', 'dropfirstlast', 'randstartend'

function [ind, siz] = sampling(sz, samplevalues, method)
if strcmp(method, 'oldmethod')
    ind = [];
    siz = [];
    for sp=samplevalues
        inds = 1:sp:sz;        
        ind = [ind inds];
        siz = [siz sp*ones(1, length(inds))];
    end
    return
end

if strcmp(method, 'randstartend')
    edges = sort(round(rand(1,2)*sz));
    sz = sz - sum(edges);
end

rest = sz;%should be -1 only -2 for the dirty 1nm fix
position = 1;
ii = 1;

while rest >= samplevalues(1)
    ind(ii) = position;
    rnd = ceil(sum(samplevalues<=rest)*rand(1));
    
    if(rnd == 0)
        rnd = 1;
    end
    siz(ii) = samplevalues(rnd);
    
    rest = rest-siz(ii);
    position = position+siz(ii);
    ii = ii+1;
end

if ~exist('ind')
    ind=[];
    siz=[];
end

if ~strcmp(method, 'pla')
    rnd_ind = randperm(size(siz,2));
    siz = siz(rnd_ind);
    ind = [1 cumsum(siz)+1];
end

if strcmp(method, 'dropfirstlast')
    % remove first and last segments
    ind = ind(2:end-1);
    siz = siz(2:end-1);
end

if strcmp(method, 'randstartend')
    ind = ind + edges(1);
end

