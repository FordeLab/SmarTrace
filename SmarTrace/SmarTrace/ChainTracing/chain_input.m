function [ X, Y ] = chain_input()
% Copyright 2016 Naghmeh Rezaei
% ******Do NOT distribute******
X=[];
Y=[];
cX=1;
cY=1;
plt_handles = {};

    function clear_plot()
        for hh=plt_handles
            try
                delete(hh{1});
            catch E
                %
            end
        end
        plt_handles = {};
    end

    function plot_now()
        clear_plot();
        %plt_handles{end+1} = plot(X, Y, '--o', 'Color', [0 0 0]);
        plt_handles{end+1} = plot(X, Y, 'o', 'Color', [0 0 0],'markersize',14);
        if length(X)>1
            [sp_val, dir_der] = mod_interparc(1e-6:2:600, X, Y, 'spline');
            plt_handles{end+1} = plot(sp_val(:,1), sp_val(:,2), '--', 'Color', [0 0 0],'linewidth',1.5);            
        end
    end

    function mouseMove(object, eventdata)
        C = get (gca, 'CurrentPoint');
        cX = C(1,1);
        cY = C(1,2);
        if length(X)<1
            return
        end
        %title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
        x=xlim;
        if cX<x(1)
            xlim(x-(x(2)-x(1))*0.01)
        end
        if cX>x(2)
            xlim(x+(x(2)-x(1))*0.01)
        end
        y=ylim;
        if cY<y(1)
            ylim(y-(y(2)-y(1))*0.01)
        end
        if cY>y(2)
            ylim(y+(y(2)-y(1))*0.01)
        end
    end

    function mouseDown(object, eventdata)
        C = get (gca, 'CurrentPoint');
        %disp(['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
        X(end+1) = C(1,1);
        Y(end+1) = C(1,2);
        plot_now();
    end

    function zoomat(z)
        %axis auto
        xlim((xlim-cX)*z+cX);
        ylim((ylim-cY)*z+cY);
        axis square
    end

    function keyPress(object, eventdata)
        %eventdata
        if strcmp(eventdata.Key, 'return') ...
            || strcmp(eventdata.Key, 'c')
            uiresume
        end
        if strcmp(eventdata.Character, 'd')
            X = X(1:end-1);
            Y = Y(1:end-1);
            plot_now();
        end
        if strcmp(eventdata.Character, '-')
            zoomat(2);
        end
        if strcmp(eventdata.Character, '=')
            zoomat(0.5);
        end
        if strcmp(eventdata.Character, 'a')
            axis auto;
        end
        if strcmp(eventdata.Character, 'z')
            %zooms in around the selected chain
            X1=min(min(X)); X2=max(max(X)); Xc = (X1+X2)/2;
            Y1=min(min(Y)); Y2=max(max(Y)); Yc = (Y1+Y2)/2;
            d = max([X2-X1 Y2-Y1])/2+10;
            axis([ Xc-d Xc+d Yc-d Yc+d]);
            axis square;
        end        
        if strcmp(eventdata.Character, 'q') ...
            || strcmp(eventdata.Key, 'escape')
            X=[];
            Y=[];
            uiresume
        end
    end

%figure(gcf);
set(gcf, 'WindowButtonMotionFcn', @mouseMove);
set(gcf, 'WindowButtonDownFcn', @mouseDown);
set(gcf, 'KeyPressFcn', @keyPress);
set(gcf, 'Pointer', 'cross')
ish = ishold();
hold on;
uiwait;
clear_plot; 
if (~ish)
    hold off;
end
set(gcf, 'Pointer', 'arrow')
set(gcf, 'WindowButtonMotionFcn', '');
set(gcf, 'WindowButtonDownFcn', '');
set(gcf, 'KeyPressFcn', '');
drawnow;

end

