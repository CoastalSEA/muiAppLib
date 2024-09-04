function [E,N,dist,Beta,Alpha] = crenulate_bays(inp)
%
%-------function help------------------------------------------------------
% NAME
%   crenulate_bays.m
% PURPOSE
%   Generate the shoreline for an equilibrium crenulate bay using the
%   method of Hsu and Evans, 1989
% USAGE
%   [E,N,dist,Beta,Alpha] = crenulate_bays()
% INPUTS
%   inp - defines intial default values to use if empty or not used
%   UI interface used to input variables
% OUTPUT
%   E, N - vectors of Eastings and Northing of the crenulate shore,
%   dist - length of control line (m)
%   Alpha - angle of wave crest to True North (deg TN) - clockwise from TN
%           i.e. 90 deg is a wave propagting to the south
%   Beta - angle between wave crest and contol line (deg)
% NOTES
%   Hsu J R C and Evans C, 1989, Parabolic bay shapes and applications. 
%   Proceedings - Institution of Civil Engineers.Part 2.Research and 
%   theory, 87, 557-570.
%   Hsu J R C, Silvester R and Xia Y M, 1989, Static equilibrium bays: new 
%   relationships. Journal of Waterway, Port, Coastal and Ocean Engineering, 
%   ASCE, 115 (3), 285-298.
%
% Written by Ivan Haigh and Ian Townend
% 29th July 200, Nov 2018, Mar 2022
%--------------------------------------------------------------------------
%
    ok=1;
    if nargin<1 || isempty(inp)
        %default data
        Ro = 5000;        %length of control line
        beta = 45;        %Angle between control line and wave crest
        Eo = 258000;      %Easting of control point
        No = 300000;      %Northing of control point
        alpha = 90;       %Angle of wave crest to TN
        p = 0;            %Spiral rotation (1=clockwise outwards and 0=counter-clockwise) 
        isplot = true;    %create figure for plot
        ishold = false;   %logical flag, true to keep any previous shoreline
    else
        Ro = inp.Ro;      %length of control line
        beta = inp.beta;  %Angle between control line and wave crest
        Eo = inp.Eo;      %Easting of control point
        No = inp.No;      %Northing of control point
        alpha = inp.alpha;%Angle of wave crest to TN
        p = inp.p;        %Spiral rotation (1=clockwise outwards and 0=counter-clockwise)
        isplot = false;   %use inp.ax instead, so plot is handled by calling function
        ishold = inp.hold;%logical flag, true to keep any previous shoreline
    end
    
    %STEP 1 - GET INPUT DATA    
    prompt = {'Length of control line (m)',...
              'Angle between control line and wave crest (1 to 89deg)',...
              'Easting of control point',...
              'Northing of control point',...
              'Angle of wave crest to TN (deg)',...
              'Spiral rotation about control point (1=clockwise and 0=anticlockwise)',...
              'Keep current line (0/1)'};  
    title = 'Crenulate Bay';
    numlines = 1;

    while ok>0
        %compute bay co-ordinates and plot
        [E,N,dist,Beta,Alpha] = bayparameters(Ro,beta,Eo,No,alpha,p);
        if ~isempty(E) 
            if isplot      
                plotbay([],E,N,Eo,No,dist,alpha,beta,p);
            elseif isfield(inp,'ax') && ~isempty(inp.ax)
                if ~ishold
                    hlines = findobj(inp.ax,'Tag','c_bay');
                    delete(hlines)
%                 else
%                     hlines = findobj(inp.ax,'-regexp','Tag','c_bay*');
                end
%                 delete(hlines)
                hold on
                plotbay(inp.ax,E,N,Eo,No,dist,alpha,beta,p);
                hold off
            end
        end
        %prompt user to change default values
        defaultvalues = {num2str(Ro),num2str(beta),num2str(Eo),num2str(No),...
                            num2str(alpha),num2str(p),num2str(ishold)};
        useInp = inputdlg(prompt,title,numlines,defaultvalues);
        if isempty(useInp)
            ok = 0;
        else   
            Ro = str2double(useInp{1});
            beta = str2double(useInp{2});
            Eo = str2double(useInp{3});
            No = str2double(useInp{4});
            alpha = str2double(useInp{5});
            p = str2double(useInp{6});
            ishold = logical(str2double(useInp{7}));
        end
    end
end
%%
function [E,N,dist,angle,Alpha] = bayparameters(Ro,beta,Eo,No,alpha,p)
    %STEP 2 - CALCULATE BAY PARAMETERS 
    if beta<1 || beta>89
        warndlg('Wave angle is out of bay. Change angle between control line and wave crest')
        E = []; N = []; dist = []; angle = []; Alpha = [];
        return;
    end
    
    npts = 4200; %number of points in coastline
    nint = 0.05; %angle increment used to compute coastline
    
    %(Values from table 2 of HSU and Evans, 1989 for given beta.
    C = lookuptable();
    % a = Ro*((0.014*beta) - (0.000094*(beta^2)));

    C0 = interp1(C(:,1), C(:,2), beta, 'linear', 'extrap');
    C1 = interp1(C(:,1), C(:,3), beta, 'linear', 'extrap');
    C2 = interp1(C(:,1), C(:,4), beta, 'linear', 'extrap');

    i = 1:1:npts;
    theta = beta + nint*i;

    R1 = Ro.*( C0 + (C1.*(beta./theta)) + (C2.*((beta./theta).^2)) );

    X = R1.*cosd(theta);
    Y = R1.*sind(theta);

    if p~=1
        X = -X;
    end
    
    %STEP 3 - COVERT TO EASTINGS AND NORTHINGS
    N = No - ( (X.*(cosd(alpha))) + (Y.*(sind(alpha))) );
    E = Eo + ( (Y.*(cosd(alpha))) - (X.*(sind(alpha))) );

    e = E(1);              %shoreline start point
    n = N(1);

    dist  = (sqrt( (Eo-e).^2 + (No-n).^2 ));  %check on Ro

    angle = atand( ((Eo-e))/((No-n)) );       %check on beta
    if Eo>=e && No>=n
        Alpha = angle-beta;
    elseif Eo>=e && No<n
        Alpha = 90+(90-abs(angle))-beta;
    elseif Eo<e && No<n
        Alpha = (180+abs(angle))-beta;
    elseif Eo<e && No>=n
        Alpha = (270+(90-abs(angle)))-beta;
    else
        Alpha = [];
    end
end
%%
function plotbay(hax,E,N,Eo,No,dist,alpha,beta,p)
    %create a wave crest line and plot results
    if isempty(hax)
        hfig = figure('Name','Crenulate-bay','Tag','PlotFig');
        hax = axes(hfig);
    end
    
    if (N(1)>=No && p==0) || (N(1)<No && p==1)
        Nw = No+sign(N(1)-No)*dist*cosd(alpha)/2;
    else
        Nw = No-sign(N(1)-No)*dist*cosd(alpha)/2;        
    end
    
    if (E(1)>=Eo && p==0) || (E(1)<Eo && p==1)
        Ew = Eo+sign(E(1)-Eo)*dist*sind(alpha)/2;
    else
        Ew = Eo-sign(E(1)-Eo)*dist*sind(alpha)/2;
    end
    
    %plot results
    plot(hax,E,N,'-r','LineWidth',1.5,'DisplayName','Shoreline','Tag','c_bay');  
    hold on
    plot(hax,[E(1),Eo],[N(1),No],'--g','LineWidth',1,'DisplayName','Control line','Tag','c_bay');
    plot(hax,[Ew,Eo],[Nw,No],'--.b','LineWidth',1,'DisplayName','Wave crest','Tag','c_bay');
    h1 = plot(hax,Eo,No,'-or','DisplayName','Control pint','Tag','c_bay');                           %control point
    set(get(get(h1,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); %exclude points from legend
    xlabel('Easting (m)');
    ylabel('Northing (m)');
    grid on
    axis equal
    legend
    if p==1, rotation = 'Clockwise'; else, rotation = 'Anticlockwise'; end
    title(sprintf('Wave angle: %g; Angle to control line: %g; %s',alpha,beta,rotation))
    hold off
end
%%
function C = lookuptable()
    %lookup matrix. Values from table 2 of HSU and Evans, 1989 for given beta.
    C = [0,0,1,0;...
     5,0.018,1.005,-0.023;...
     10,0.036,1.011,-0.047;...
     15,0.050,0.998,-0.049;...
     20,0.055,1.029,-0.088;...
     25,0.054,1.0830,-0.142;...
     30,0.045,1.146,-0.194;...
     35,0.029,1.220,-0.253;...
     40,0,1.326,-0.332;...
     45,-0.039,1.446,-0.412;...
     50,-0.088,1.588,-0.507;...
     55,-0.151,1.756,-0.611;...
     60,-0.227,1.930,-0.706;...
     65,-0.315,2.113,-0.800;...
     70,-0.409,2.284,-0.873;...
     75,-0.505,2.422,-0.909;...
     80,-0.600,2.520,-0.906]; 
end