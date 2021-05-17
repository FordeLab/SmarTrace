function varargout = SmarTrace(varargin)

% copyright 2010-2016 Guillaume Lamour, Naghmeh Rezaei
% *******Do NOT distribute*******
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Edit by Mathew Schneider:
%Line 1005 - properly saves the .curv structure in handle
%Edit by Alaa Al-Shaer
% Line 630 - adds chain numbers 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Edit the above text to modify the response to help SmarTrace


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SmarTrace_OpeningFcn, ...
    'gui_OutputFcn',  @SmarTrace_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end



% --- Executes just before SmarTrace is made visible.
function SmarTrace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SmarTrace (see VARARGIN)

ifield = 'chain';
ivalue = {[];[];[]};
global directData;
directData = struct(ifield,ivalue);

% setup part of the closing confirmation pop-up window (see last function)
set(handles.figure1,'CloseRequestFcn',@closeGUI);

% masks the axis of the image
axis off
% Choose default command line output for SmarTrace
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SmarTrace wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SmarTrace_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in btn_browse_input.
function btn_browse_input_Callback(hObject, eventdata, handles)

% grab the file containing the height maps
% input_filename = uigetfile('*.txt')
[input_filename,PathName] = uigetfile('*','Select the height map');
%fid = fopen(input_filename)
%[input_filename,rt,n] = fopen(fid)

load_imagefile(input_filename, PathName, hObject, handles);

% --- Executes on button press in btn_browse_input.
function load_imagefile(input_filename, PathName, h, handles)
    if input_filename == 0
        return
    end
    
    set(handles.ans_filename, 'String', input_filename);
    handles.ans_filepath = PathName;

    [~, ~, ext] = fileparts(input_filename);
    handles.ext = ext;

    % switch directory
    cd(PathName)
    % read the heights map
    switch ext
        case '.txt'
            A = dlmread(input_filename, '\t', 1, 1);
        case {'.bmp','.png','.jpg','.ppm','.tif'}
            A = imread(input_filename);
            if size(A, 3) > 1
                A = rgb2gray(A(:,:,1:3));
            end
        otherwise
            msgbox(sprintf('This software accepts only files with following extensions: .bmp, .png, .jpg, .ppm, .tif, and .txt.  Please convert your image in one of these associated formats.'),...
                'Wrong Image File Format','Warn')
            return
    end
    
    % get resolution and store it
    image_resol1 = size(A,1);
    handles.image_resol = image_resol1;
    set(handles.txt_resol1,'String',image_resol1);
    image_resol2 = size(A,2);
    set(handles.txt_resol2,'String',image_resol2);
    Z = [image_resol1 image_resol2];
    tolerance = max(Z)/100;
%     if abs(image_resol1 - image_resol2) > tolerance
%         msgbox(sprintf('Please crop your height map so that both image sides share the same number of pixels (or at least that they do not differ by more than 1%%), then reload it.'),...
%             'Tolerance threshold exceeded','Warn')
%         %return
%     else
        % normalize heights to make the image readable
        C = round((A-min(min(A)))./(max(max(A))-min(min(A)))*255-min(min(A)));
        
        switch ext
            case '.txt'
                image(C);
            otherwise
                image(A);
        end
        colormap(bone(255));
        axis off
        axis equal % NRtest
        % Store filename
        handles.input_filename = input_filename;
        handles.A = A;
        handles.C = C;
%     end
    guidata(h,handles);


% --------------------------------------------------------------------
% Select fibril
function btn_select_Callback(h, eventdata, handles, directData)
getpoints(h,handles);


function d = angdiff(th1, th2)
d = th1 - th2;
d = mod(d+pi, 2*pi) - pi;

% --------------------------------------------------------------------
% Gets points on the fibril from user input, traces the
% fibril that corresponds to these points, and fits a spline
% to the resulting points.
function getpoints(h,handles,directData)
try
prf = 0;

ext = handles.ext;

text1 = 'Add next';
set(handles.txt_confirm_add_chain,'String',text1,...
    'foregroundcolor',[1 0.5 0],'FontWeight','bold');

C = handles.C;
A = handles.A;

colormap(bone(255));
axis off;
hold on;


% if the chain already traced, it delets the exiting trace and starts over. 
chain_n = str2num(get(handles.txt_no,'String'));
if (isfield(handles, 'chain_plots'))
    %display({'NRtest CHAIN', size(handles.chain_plots,1), chain_n})
    if size(handles.chain_plots,1)>=chain_n
        for hhh=handles.chain_plots(chain_n,:)
            try
                delete(hhh{1});
            catch E
                %
            end
        end
    end
else

    %initializes the list of plotted chains
    handles.chain_plots=cell(0,0); 
end

plot_nn = 1;


% Collect points on a fibril
% gathers an unlimited number of points until the return key is pressed
[X, Y] = chain_input;

% X = [9.1309   22.1152   31.9581   34.6806   36.6702   42.9529   51.5393];
% 
% 
% Y = [15.1073   21.1806   30.0812   34.4791   41.1806   45.3691   45.6832];
   
if length(X)<2
   disp('not enough points select') 
   return
end
    
if prf
    profile clear 
    profile on
end


set(gcf, 'Pointer', 'watch');

%zooms in around the selected chain
X1=min(min(X)); X2=max(max(X)); Xc = (X1+X2)/2;
Y1=min(min(Y)); Y2=max(max(Y)); Yc = (Y1+Y2)/2;
d = max([X2-X1 Y2-Y1])/2+10;
axis([ Xc-d Xc+d Yc-d Yc+d]);
axis square;
hold on;

drawnow;

%draw scale bar
image_width = str2double(get(handles.txt_input_width,'String'));
nmperpx=image_width;
% x1=max(X);
% y1=max(Y)+10;
% scale=20.0/nmperpx;
% plt = line([x1-scale x1], [y1 y1], 'LineWidth', 10, 'Color', [1 1 1]);
% handles.chain_plots{chain_n,plot_nn} = plt;
% plot_nn = plot_nn+1;

% Show user defined points
plt = plot(X, Y, 'ko', 'markersize', 10);
handles.chain_plots{chain_n,plot_nn} = plt;
plot_nn = plot_nn+1;

%b-spline fit to the user defined points 
%returns spline values (sp_val) and its derivatives (dir_der)
%modified interparc fuction fit returns equially distanced points
%note dir_der is not unit length
[sp_val, dir_der] = mod_interparc(1e-6:600,X,Y,'spline'); 
Xsp = sp_val(:,1);
Ysp = sp_val(:,2);


%initial b_spline to the user defiend points in nm:
l2 = length(Xsp)*nmperpx*1.2;
[sp_val1, dir_der1,xxx, kappa1] = mod_interparc2(1e-6:l2,X*nmperpx,Y*nmperpx,'spline');

Xsp1 = sp_val1(:,1);
Ysp1 = sp_val1(:,2);

handles.Xsp1 = Xsp1;
handles.Ysp1 = Ysp1;
handles.dir_der1 = dir_der1;
handles.kappa1 = kappa1;

% Remove background
se = strel('disk', 12); %default 12
im2 = imtophat(A, se);

% median filter to remove speckles, use odd sizes to prevent shifting
im1 = medfilt2(im2, [3 3]);

[hh ww] = size(im1);
imcenter = im1(round(hh*0.1):round(hh*0.9),round(ww*0.1):round(ww*0.9));
THRESH = double(prctile(reshape(imcenter,1,[]), 80));
MAX = double(prctile(reshape(imcenter,1,[]), 99.9));
% display(sprintf('Using Threshold=%f Max=%f', THRESH, MAX));

% plt = imshow(im1);
% handles.chain_plots{chain_n,plot_nn} = plt;
% plot_nn = plot_nn+1;
% drawnow;

% Vectors for storing the traced points
% NRtest count = 0;

for count=1:1
    
XT = Xsp;
YT = Ysp;

% searh range
steps = 0.25; %default 0.2
trange = -6:steps:6; %default 3
prange = -2:steps:2; %default -2

tr=repmat(trange,length(prange),1);
pr=repmat(prange',1,length(trange));


lwcv = zeros(0,0,0); % value (score) at a specific length, width & centre compare to where the initial spline is.

for ll=1:length(Xsp)
    %     plt = plot(Xc, Yc, '*r');
    %     handles.chain_plots{chain_n,plot_nn} = plt;
    %     plot_nn = plot_nn+1;
    
    % tangent direction
    dx = dir_der(ll,1); % dir-der from mod_interparc
    dy = dir_der(ll,2);
    dl = sqrt(dx^2+dy^2);
    dx = dx / dl; % normal tangent vector  
    dy = dy / dl;
    
    % can use meshgrid instead
    Px = Xsp(ll) + tr * dy + pr * dx; % x corrdinates of search area (grid) 
    Py = Ysp(ll) + tr * -dx + pr * dy;
    
    % remove out-of-bound points
    Px(Px<1)=1;
    Py(Py<1)=1;
    Px(Px>size(im1,1)-1)=size(im1,1)-1; % repeats border px for out of bound coords
    Py(Py>size(im1,2)-1)=size(im1,2)-1;
    
    %     if ll == 5
    %         plt = plot(Px, Py, 'sg');
    %         handles.chain_plots{chain_n,plot_nn} = plt;
    %         plot_nn = plot_nn+1;
    %     end
    
    pts = [];
    %         pts = get_subpixel(im1, [Px ; Py]'-1, 'cubic')';
    
    % % TODO: remove this part of the code, sub pixel calculations
    Pxflat = reshape(Px,1,[]); %grid x points into a list
    Pyflat = reshape(Py,1,[]);
    % NRtest 
    pts = pkget_subpixel(im1, [Pxflat ; Pyflat]', 'cubic')'; 
%     pts = pkget_subpixel(im1, [Pxflat ; Pyflat]', 'nearest')'; 
    pts = reshape(pts,length(prange),length(trange)); % intensities at the grid points
    
    % % direct pixel calculations
    %inds = sub2ind(size(im1), round(Py), round(Px));
    %pts = im1(inds);
    %pts = reshape(pts,length(prange),length(trange));
    
        
    if ll == 20
        wcv = centroid06(pts, 0, THRESH, MAX);
    else
        wcv = centroid06(pts, 0, THRESH, MAX);
    end
        
    
    if (max(size(lwcv))==0)
        lwcv=zeros(length(Xsp),size(wcv,1),size(wcv,2));
    end
    lwcv(ll,:,:,:) = wcv;
        
    %             ff = gcf;
    %             figure(2);
    %             plot(mat2gray(pts),'b-');
    %             hold on;
    % %             plot(mat2gray(pp),'g-');
    %             plot(ii,0.5,'r*');
    %             hold off;
    %             pause(2);
    %             figure(ff);
    
    %             v = pp/4+0.5;
    %             plt = plot(Px,Py,'*','Color', [v ; v ; v]);
    %             handles.chain_plots{chain_n,plot_nn} = plt;
    %             plot_nn = plot_nn+1;
    
    %     t = trange(ii);
    %
    %     if 0
    %         % HACK: give more weight to elements that are near the user
    %         % selected chain
    %         dist = sqrt(Px-Xsp(ll)).^2+(Py-Ysp(ll)).^2;
    %         pts = pts .* gaussmf(dist, [15 0]);
    %
    %         % HACK: give more weight to elements that have the same direction
    %         % as the selected nodes
    %         p_dirs = atan2(Py-Yc,Px-Xc);
    %         p_dir = atan2(dir_der(ll,2),dir_der(ll,1));
    %         angle_diffs = abs(angdiff(p_dirs,p_dir));
    %         pts = pts .* gaussmf(angle_diffs, [3 0]);
    %
    %         [maxv maxi] = max(pts);
    %         t = trange(maxi);
    %     end
    %
    %     XT(ll) = Xsp(ll) + dy * t;
    %     YT(ll) = Ysp(ll) + -dx * t;
    %
    %     %         plt = plot(PcenterX, PcenterY, '*g');
    %     %         handles.chain_plots{chain_n,plot_nn} = plt;
    %     %         plot_nn = plot_nn+1;
    %
    %     XJ = [XT(ll)-dy*steps*ww/2 XT(ll)+dy*steps*ww/2];
    %     YJ = [YT(ll)+dx*steps*ww/2 YT(ll)-dx*steps*ww/2];
    %     plt = plot(XJ, YJ, 'om','markersize', 3);
    %     handles.chain_plots{chain_n,plot_nn} = plt;
    %     plot_nn = plot_nn+1;
    
    %     plt = plot(XT(ll), YT(ll), '*', 'markersize', 10, 'color', [max(0,min((vv-0.3)*3,1)) 0 0]);
    %     handles.chain_plots{chain_n,plot_nn} = plt;
    %     plot_nn = plot_nn+1;
    
end

% for ll=1:length(Xsp)
%     wcv=lwcv(ll,:,:);
%     [v,idx] = max(wcv(:));
%     [a,w,c] = ind2sub(size(wcv),idx);
%     ii=c;
%     ww=w;
%
%     dx = dir_der(ll,1);
%     dy = dir_der(ll,2);
%     dl = sqrt(dx^2+dy^2);
%     dx = dx / dl;
%     dy = dy / dl;
%
%     t = trange(ii);
%     XT(ll) = Xsp(ll) + dy * t;
%     YT(ll) = Ysp(ll) + -dx * t;
%
%     %         plt = plot(PcenterX, PcenterY, '*g');
%     %         handles.chain_plots{chain_n,plot_nn} = plt;
%     %         plot_nn = plot_nn+1;
%
%     XJ = [XT(ll)-dy*steps*ww/2 XT(ll)+dy*steps*ww/2];
%     YJ = [YT(ll)+dx*steps*ww/2 YT(ll)-dx*steps*ww/2];
%     plt = plot(XJ, YJ, 'om','markersize', 3);
%     handles.chain_plots{chain_n,plot_nn} = plt;
%     plot_nn = plot_nn+1;
% end

lwcv(lwcv<0.01)=0.01;

% chain points
chpX = zeros(length(Xsp),2);
chpY = zeros(length(Xsp),2);

% find widths with max values (scores) along chain
for ll=1:length(Xsp)
    wcv = squeeze(lwcv(ll,:,:));
    
    %findpeaksn finds local maxima and ignores nearby peaks
    [cw,cc]=ind2sub(size(wcv), find(findpeaksn(wcv,true(1, ndims(wcv)), 0.1, 5)));
    if length(cw) == 0
        warning('no peaks at length=%d', ll);
        wl(ll) = 3/steps; % NRtest wl(ll-1);
        continue;
    end
    wl(ll) = mean(cw);
    
    %     idx=find(lwcv(ll,:,:)>0.2);
    %     [l,w,c,v]=ind2sub(size(lwcv), idx);
    %     v=v/sum(v);
end
% smooth widths, moving average
norm_filter = fspecial('gauss', [1 10], 3);
wl = imfilter(wl, norm_filter, 'replicate');
for ll=1:length(Xsp)
    % TODO:: refactor this for faster speed
    dx = dir_der(ll,1);
    dy = dir_der(ll,2);
    dl = sqrt(dx^2+dy^2);
    dx = dx / dl;
    dy = dy / dl;
    
    % give more weight to those near the smooth width at that point of the
    % chain
    wcv = squeeze(lwcv(ll,:,:));
    W = repmat(normpdf(1:size(wcv,1), wl(ll), size(wcv,1)*0.1)', 1, size(wcv,2));
    W = W/max(max(W));
    wcv = W.*wcv;
    
    % more weight for points in the same direction as the last point
    if ll>2
        tdx = XT(ll-1) - XT(ll-2);
        tdy = YT(ll-1) - YT(ll-2);
        tdl = sqrt(tdx.^2+tdy.^2);
        tdx = tdx / tdl;
        tdy = tdy / tdl;
        
        xn = XT(ll-1) + tdx;
        yn = YT(ll-1) + tdy;
        
        Wc=[];
        for cc=1:size(wcv,2)
            c = trange(cc);
            
            x2 = Xsp(ll) + dy * c;
            y2 = Ysp(ll) + -dx * c;
            l2 = (x2-xn).^2+(y2-yn).^2;
            Wc(cc) = 5/(5+l2);
        end
        W=repmat(Wc, size(wcv,1), 1);
        wcv = W.*wcv;
    end
    
    [cw, cc] = ind2sub(size(wcv), find(findpeaksn(wcv,true(1, ndims(wcv)))));
    if length(cw) == 0
        error('No peaks found')
    end
    cv=[];
    for ii=1:length(cw)
        cv(end+1)=wcv(cw(ii),cc(ii));
    end
    %[cw';cc';cv]
    [a,k]=max(cv);
    ii=cc(k);
    ww=cw(k);
    vv=cv(k);
    
    
    %     [a,k]=max(wcv(:));
    %     [cw,cc,cv] = ind2sub(size(wcv),k);
    %     ww=cw;
    %     ii=cc;
    
    t = trange(ii);
    
    XT(ll) = Xsp(ll) + dy * t;
    YT(ll) = Ysp(ll) + -dx * t;
    
    XJ = [XT(ll)-dy*steps*ww/2 XT(ll)+dy*steps*ww/2];
    YJ = [YT(ll)+dx*steps*ww/2 YT(ll)-dx*steps*ww/2];
    chpX(ll,:) = XJ;
    chpY(ll,:) = YJ;
        
%     cm = hot;
%     rgb = ind2rgb(round(vv * 255), cm);
%     if ll>1
%         plt = plot(XT(ll-1:ll), YT(ll-1:ll), '-', 'LineWidth', 2, 'Color', squeeze(rgb));
%         handles.chain_plots{chain_n,plot_nn} = plt;
%         plot_nn = plot_nn+1;
%         drawnow;
%     end
end

% NRtest
cm = lines;
rgb = squeeze(ind2rgb(round(vv * 255), cm));
% if plot_nn<5
%     plt = plot(XT, YT, 'r--', 'LineWidth', 1,'color',rgb);
%     handles.chain_plots{chain_n,plot_nn} = plt;
%     plot_nn = plot_nn+1;
% else
%     plt = plot(XT, YT, 'r--', 'LineWidth', 2,'color',rgb);
%     handles.chain_plots{chain_n,plot_nn} = plt;
%     plot_nn = plot_nn+1;
%     drawnow
% end
Xsp = XT;
Ysp = YT;

%pause(0.5);

end

plt = plot(chpX(:,1), chpY(:,1), '--b', 'LineWidth', 1.5);
handles.chain_plots{chain_n,plot_nn} = plt;
plot_nn = plot_nn+1;
plt = plot(chpX(:,2), chpY(:,2), '--b', 'LineWidth', 1.5);
handles.chain_plots{chain_n,plot_nn} = plt;
plot_nn = plot_nn+1;

set(gcf, 'Pointer', 'arrow')


XT=XT'; %the final traced points on the chain
YT=YT';

%tpoints = [1:length(XT)];
%save('trialfile','XT','YT','tpoints');

% plt = plot(XT, YT, 'r.');
% handles.chain_plots{chain_n,plot_nn} = plt;
% plot_nn = plot_nn+1;

handles.X = X;
handles.Y = Y;
handles.XT = XT;
handles.YT = YT;

% ---------------------------------------------------------------
% --- Compile spline fitting and image data (--> normalize lengths vs. pixels)
%---- Executes on button press in btn_compile_and_check.
image_width = str2double(get(handles.txt_input_width,'String'));
handles.image_width = image_width;
image_resol = handles.image_resol;

output_base = get(handles.txt_output_path,'String');
output_no = get(handles.txt_no,'String');
file_name=[output_base, output_no];

nmperpx = image_width;

%nmperpx = image_width*1000/image_resol % lengths are converted from microns to nanometers (*1000)

DOTS = [XT ; YT];

% DOTS = DOTS(:,find(sum(diff(DOTS').^2,2)~=0));
l2 = length(Xsp)*nmperpx*1.2;
% [sp_val, dir_der] = mod_interparc(1e-6:l2,DOTS(1,:)*nmperpx,DOTS(2,:)*nmperpx,'spline'); %note dir_der is not unit length

%return curvature, derivatives and spline
[sp_val, dir_der, xxx, kappa] = mod_interparc2(1e-6:l2,DOTS(1,:)*nmperpx,DOTS(2,:)*nmperpx,'spline');

x_norm = sp_val(:,1);
y_norm = sp_val(:,2);

Xsp1 = handles.Xsp1;
Ysp1 = handles.Ysp1;

plt = plot(Xsp1/nmperpx, Ysp1/nmperpx, 'k--','linewidth',1);
handles.chain_plots{chain_n,plot_nn} = plt;
plot_nn = plot_nn+1;

plt = plot(x_norm/nmperpx, y_norm/nmperpx, 'r-','linewidth',1.5);
%add chain number
text(x_norm(1)/nmperpx, y_norm(1)/nmperpx, num2str(chain_n),'HorizontalAlignment', 'right', 'Color', 'g');
handles.chain_plots{chain_n,plot_nn} = plt;
plot_nn = plot_nn+1;



handles.dir_der = dir_der;
handles.kappa = kappa;

% plt = text(max(x_norm)+5, max(y_norm)+5, num2str(chain_n), 'Color', 'black', 'BackgroundColor', 'white', 'Margin', 0.1, 'fontsize', 14);
% handles.chain_plots{chain_n,plot_nn} = plt;
% plot_nn = plot_nn+1;

hold off;

bspline_norm_coord(1,:) = x_norm;
bspline_norm_coord(2,:) = y_norm;
handles.bspline_norm_coord = bspline_norm_coord;


handles.file_name = file_name;
handles.output_base = output_base;
handles.output_no = output_no;

handles.x_norm = x_norm; % used in interval/check increment function
handles.y_norm = y_norm;

handles.nmperpx = nmperpx;

arclen = arclength(x_norm,y_norm)

% workingXT = strcat('XT_',num2str(chain_n));
% workingYT = strcat('YT_',num2str(chain_n));
% S.(workingXT) = XT;
% save('trialfile.mat','-struct','S','-append');
% S.(workingYT) = YT;
% save('trialfile.mat','-struct','S','-append');
% tpoints = [0:length(XT)-1];
% workingT = strcat('tpoints_',num2str(chain_n));
% S.(workingT) = tpoints;
% save('trialfile.mat','-struct','S','-append');

global directData;

field = 'chain';
tpoints = [0:length(XT)-1];
value = {XT;YT;tpoints};
tempStruct = struct(field,value);
directData = [directData,tempStruct];



% ---------------------------------------------------------------
% Check knots spline increments over the fibril
%INTERVAL FUNCTION
x_norm = handles.x_norm;
y_norm = handles.y_norm;

for i=1:length(x_norm)-1
    dist = ( x_norm(i) - x_norm(i+1) ).^2 + ( y_norm(i) - y_norm(i+1) ).^2;
    intervals(i) = sqrt(dist);
end

% average_interval between each knot over all the spline
mean_itv = mean (intervals(1,:));
mean_itv = round (mean_itv);
% standard deviation
sd_itv = std(intervals(1,:));
sd_itv = round (sd_itv);

%-------- Send the values in the box of the GUI dedicated to it
% set(handles.ans_mean_itv,'String',mean_itv);
% set(handles.ans_sd_itv,'String',sd_itv);

handles.mean_itv = mean_itv;
handles.sd_itv = sd_itv;
handles.intervals = intervals;

if prf
    profile viewer
    %profile off
end

guidata(h,handles);
catch E
guidata(h,handles);
rethrow(E)
end


function txt_fitting_param_Callback(hObject, eventdata, handles)
%store the contents of input1_editText as a string. if the string
%is not a number then input will be empty
input = str2num(get(hObject,'String'));
%checks to see if input is empty. if so, default input1_editText to twenty
if (isempty(input))
    set(hObject,'String','15')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_fitting_param_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
% --- user input - image side size in microns
function txt_input_width_Callback(hObject, eventdata, handles)
% hObject    handle to txt_input_width (see GCBO)

function txt_input_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Setup part of filename by user input
function txt_output_path_Callback(hObject, eventdata, handles)
% hObject    handle to txt_output_path (see GCBO)

% --- Executes during object creation, after setting all properties.
function txt_output_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Setup fibril # to be saved by user input
function txt_no_Callback(hObject, eventdata, handles)
% hObject    handle to txt_no (see GCBO)

% --- Executes during object creation, after setting all properties.
function txt_no_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%---- Executes on button press in btn_combocalc.
%-------Executes all three procedures at the same time
function btn_combocalc_Callback(h, eventdata, handles)
% hObject    handle to btn_combocalc (see GCBO)
combo_calc(h, handles);

function combo_calc(h, handles)


%TANTAN-COREL
%--------------------------------------------------------------------
%-------Calculate the decay of tangent-tangent correlations over the fibril
x_norm = handles.x_norm;
y_norm = handles.y_norm;
nmperpx = handles.nmperpx;

arclen = arclength(x_norm,y_norm);

seps_mean = [];
cosines_mean = [];
cosines_std = [];

R2s_mean = [];
R2s_std = [];

cosines = [];
seps = [];
R2s = [];

angles2pi = [];
seps_angles = [];

thetas = [];


separations = 1:1:max(arclen);
for sep = separations
    
    %ll = 0:sep:max(arclen);
    %[x_points, y_points, xder,yder, idx_points] = point_at_length(handles,ll);
    
    x_points = handles.x_norm(1:sep:end);
    y_points = handles.y_norm(1:sep:end);
    der = normr(handles.dir_der(1:sep:end,:));
    
    %correlation (cos(tetha)) calculations
    %der = [xder, yder];
    
    raw_theta = (atan2(der(:,2),der(:,1)));
    
    theta = abs(diff(raw_theta));
    thetas = [thetas;theta];
    
    angle = abs(mod(diff(raw_theta)+pi,2*pi)-pi);
    angles2pi = [angles2pi;angle];
    seps_angles = [seps_angles;ones(length(angle),1)*sep];
    
    cosine = dot(der(1:end-1,:),der(2:end,:),2);
    cosines = [cosines;cosine];
    seps = [seps;ones(length(cosine),1)*sep];
    cosines_mean = [cosines_mean;mean(cosine)];
    cosines_std = [cosines_std;std(cosine)];
    seps_mean = [seps_mean,sep];
    
    %end-to-end (R^2) calculations
    rr = [x_points,y_points];
    R2 = sum(diff(rr).^2,2);
    R2s = [R2s;R2];
    
    R2s_mean = [R2s_mean;mean(R2)];
    R2s_std = [R2s_std;std(R2)];
end

handles.angles2pi = angles2pi;
handles.cosines = cosines;
handles.contours = seps;
handles.thetas = thetas;

corel_single = [cosines,seps];
handles.corel_single = corel_single;

angles2pi_single = [angles2pi,seps_angles];
handles.angles2pi_single = angles2pi_single;

thetas_single = [thetas,seps_angles];
handles.thetas_single = thetas_single;

handles.end2end = R2s;

end2cont_single = [R2s,seps];
handles.end2cont_single = end2cont_single;

goodcontour = arclen;
handles.goodcontour = goodcontour ;

handles.total_end2end = sqrt (( y_norm(end) - y_norm(1) ).^2 +...
    (( x_norm(end) - x_norm(1) ).^2 ));


%--------------------------------------------------------------------
%BELOW is the code that fill in the structure with all the elements, pertaining
% to one fibril, calculated above in this very same function

% IMPORTANT to keep this
output_no = handles.output_no;
n=str2num(output_no);

handles.mat_length_struct(n).contourL = handles.goodcontour;
handles.mat_intervals_struct(n).meanitv_fib = handles.mean_itv;
handles.mat_sd_struct(n).meansd_fib = handles.sd_itv;

handles.bsplines_struct(n).bspline = handles.bspline_norm_coord;

handles.correlations_struct(n).corelfib = handles.corel_single;
handles.wormlike_struct(n).wormfib = handles.end2cont_single;

handles.length_conformation(n, 1) = handles.goodcontour;
handles.length_conformation(n, 2) = handles.total_end2end;

handles.points(n).X = handles.X;
handles.points(n).Y = handles.Y;
handles.points(n).XT = handles.XT;
handles.points(n).YT = handles.YT;
handles.points(n).filename = get(handles.ans_filename, 'String');
handles.points(n).filepath = handles.ans_filepath;

handles.der(n).Xder = handles.dir_der(:,1);
handles.der(n).Yder = handles.dir_der(:,2);

handles.curv(n).curvature = handles.kappa;

handles.angles2pi_struc(n).angles2pi = handles.angles2pi_single;
handles.thetas_struc(n).thetas = handles.thetas_single;

handles.info_struc(n).image_width = handles.image_width;
handles.info_struc(n).image_resol = handles.image_resol;

handles.userspline(n).Xsp1 = handles.Xsp1;
handles.userspline(n).Ysp1 = handles.Ysp1;
handles.userspline(n).dir_der1 = handles.dir_der1;
handles.userspline(n).kappa1 = handles.kappa1;


% Increment the filename
% n=str2num(output_no);
% set(handles.txt_no, 'String', num2str(n+1));

% Set the string and color of text to confirm chain has been added
text1 = 'Chain added!';
set(handles.txt_confirm_add_chain,'String',text1,...
    'foregroundcolor',[0 0.5 0],'FontWeight','bold');

guidata(h, handles);


% --- Executes on button press in loadmat.
function loadmat_Callback(hObject, eventdata, handles)

struc = uiimport

handles.mat_length_struct = struc.contour_lengths;
handles.mat_intervals_struct = struc.intervals_means;
handles.mat_sd_struct = struc.intervals_sd;

handles.bsplines_struct = struc.bsplines_norm;

handles.correlations_struct = struc.mat_correlations;
handles.wormlike_struct = struc.mat_wormlike;

handles.length_conformation = struc.length_conformation;
handles.points = struc.points;
handles.der = struc.der;
handles.angles2pi_struc = struc.angles2pi;
handles.thetas_struc = struc.thetas;
handles.curv = struc.curv;
handles.info_struc = struc.info;

handles.userspline = struc.userspline;
s = size(handles.mat_length_struct,2);
set(handles.txt_no, 'String', num2str(s+1));

guidata(hObject,handles);




% --- Executes on button press in increment_fibril_nb.
function increment_fibril_nb_Callback(hObject, eventdata, handles)
n = str2num(get(handles.txt_no,'String'));
set(handles.txt_no, 'String', num2str(n+1));
% NRtest: in future zoom to the chain here...
axis auto
axis equal

% X1=min(min(X)); X2=max(max(X)); Xc = (X1+X2)/2;
% Y1=min(min(Y)); Y2=max(max(Y)); Yc = (Y1+Y2)/2;
% d = max([X2-X1 Y2-Y1])/2+10;
% axis([ Xc-d Xc+d Yc-d Yc+d]);
% plot1 = plot(X, Y, 'bo');
% axis square;


guidata(hObject, handles);

% --- Executes on button press in pushminus.
function pushminus_Callback(hObject, eventdata, handles)
n = str2num(get(handles.txt_no,'String'));
set(handles.txt_no, 'String', num2str(n-1));
% NRtest: in future zoom to the chain here...
axis auto

guidata(hObject, handles);



% --- Executes on button press in btn_savedata.
function btn_savedata_Callback(h, eventdata, handles)

global directData;
% Important Stuff
contour_lengths = handles.mat_length_struct;
intervals_means = handles.mat_intervals_struct;
intervals_sd = handles.mat_sd_struct;

bsplines_norm = handles.bsplines_struct;

mat_correlations = handles.correlations_struct;
mat_wormlike = handles.wormlike_struct;

length_conformation = handles.length_conformation;
points = handles.points;
der = handles.der;
angles2pi = handles.angles2pi_struc;
thetas = handles.thetas_struc;

curv = handles.curv;

info = handles.info_struc;

userspline = handles.userspline;
save('dataForVectorAnalysis.mat','directData');
%output_base = handles.output_base;
output_base = get(handles.txt_output_path,'String');
output_sup = '_SmTr';
updated_filename=[output_base, output_sup];

%newfilename = fullfile(pathname, updated_filename);
save(updated_filename,'contour_lengths','intervals_means',...
    'intervals_sd','bsplines_norm',...
    'mat_correlations','mat_wormlike','length_conformation','points','angles2pi','thetas', 'info', 'der','curv','userspline');

% message confirmation
msgbox(sprintf('data successfully saved. To analyzed first run SmTranalysis1.m and then run SmTranalysis2.m.'),...
    'Confirmation Note','Help')

guidata(h, handles);



function closeGUI(src,evnt)
%this function is called when the user attempts to close the GUI window

%this command brings up the close window dialog
selection = questdlg('Do you want to close SmarTrace? Note that you can save a .mat file and reload it later, using the ''Reload Data'' button.',...
    'Close Request Function',...
    'Yes','No','No');

%if user clicks yes, the GUI window is closed, otherwise
%nothing happens
switch selection,
    case 'Yes',
        delete(gcf)
    case 'No'
        return
end


% --- Executes during object creation, after setting all properties.
function txt_nb_knots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nb_knots (see GCBO)


% --- Executes on button press in pushbut_aboutEasy.
function pushbut_aboutEasy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbut_aboutEasy (see GCBO)
msgbox(sprintf('***Do NOT distribute***\n \n 2016: Naghmeh Rezaei \n Simon Fraser University \n nra18@sfu.ca\n ----------------- \n 2010:Guillaume Lamour '),...
    'About Easyworm','Help')


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
msgbox(sprintf('contact nra18@sfu.ca.'),...
    'Help','Help')


% --- Executes on button press in btn_saveimage.
function btn_saveimage_Callback(hObject, eventdata, handles)
% hObject    handle to btn_saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
F = getframe(handles.axes4);
[im,Map] = frame2im(F);
fname = [ 'SmarTrace_' get(handles.ans_filename, 'String')];
imwrite(im, fname)


function [x,y,xder,yder,idx] = point_at_length(handles, ll)

x_norm = handles.x_norm;
y_norm = handles.y_norm;
dir_der = handles.dir_der;

if ~isfield(handles, 'length_on_chain')
    
    x_diff = diff(x_norm);
    y_diff = diff(y_norm);
    
    r_diff = sqrt(x_diff.^2 + y_diff.^2);
    
    handles.length_on_chain = cumsum(r_diff);
end
length_on_chain = handles.length_on_chain;

idx = [];
for ll2=ll
    idx(end+1) = find(length_on_chain>ll2,1);
end
x = x_norm(idx);
y = y_norm(idx);
xder = dir_der(idx,1);
yder = dir_der(idx,2);
