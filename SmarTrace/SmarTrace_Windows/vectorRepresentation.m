function vectorRepresentation(filename)
    
    numberofchains = 10;
    load(filename);
    dataContainer = struct('r2values',[],'cosvalues',[],'kappavalues',[]);
    for i = 2:numberofchains;
        XT = directData(1,i).chain;
        YT = directData(2,i).chain;
        tpoints = directData(3,i).chain;
        load(filename);
        [xout,tout] = prepareCurveData(XT,tpoints);
        [yout,tout] = prepareCurveData(YT,tpoints);
        errorvector = [];
        for c = 8:25;
            xparams = polyfit(tout,xout,c);
            yparams = polyfit(tout,yout,c);
            totalerror = 0;
            for n = 1:length(xout);
                xtemperror = (polyval(xparams,n) - xout(n))^2;
                ytemperror = (polyval(yparams,n) - yout(n))^2;
                tottemperror = xtemperror + ytemperror;
                totalerror = totalerror + tottemperror;            
            end
            errorvector = [errorvector,totalerror];
        end
        [M,i] = min(errorvector);
        xparams = polyfit(tout,xout,(i+7));
        yparams = polyfit(tout,yout,(i+7));
    
        syms t;
        xequation = (xparams(length(xparams))) + (xparams(length(xparams)-1))*t;
        yequation = (yparams(length(yparams))) + (yparams(length(yparams)-1))*t;
        for n = (length(xparams)-2):-1:1;
            xequation = xequation + (xparams(n))*t^(length(xparams)-n);
            yequation = yequation + (yparams(n))*t^(length(xparams)-n);
        end
    
        dif1xequation = diff(xequation);
        dif2xequation = diff(dif1xequation);
        dif1yequation = diff(yequation);
        dif2yequation = diff(dif1yequation);
    
        meanendsqdist = sqrt(((xequation - vpa(subs(xequation,t,1)))^2) + ((yequation - vpa(subs(yequation,t,1)))^2));
        tangentcorrelation = (dif1xequation*vpa(subs(dif1xequation,t,1)))+(dif1yequation*vpa(subs(dif1yequation,t,1)));
        curvature = abs((dif1xequation*dif2yequation)-(dif1yequation*dif2xequation))/(((((dif1xequation)^2)+(dif1yequation)^2))^(3/2));
        r2array = [];
        cosarray = [];
        kappaarray = [];
        for j = 1:0.01:(length(tpoints));
            r2array = [r2array,subs(meanendsqdist,t,j)];
            r2tempstruct = struct('r2values',r2array);
            cosarray = [cosarray,subs(tangentcorrelation,t,j)];
            costempstruct = struct('cosvalues',cosarray);
            kappaarray = [kappaarray,subs(curvature,t,j)];
            kappatempstruct = struct('kappavalues',kappaarray);
        end
        dataContainer = [datacontainer,r2tempstruct];
        dataContainer = [datacontainer,costempstruct];
        dataContainer = [datacontainer,kappatempstruct];
    end
    
    
        
            
    end
    

        
    
    