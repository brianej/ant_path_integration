function img = ImgGrabber(x,y,z,th,X,Y,Z,colp,hfov,resolution)
% Grab an image from world model
% x,y,z : camera coordinates in metres
% th : heading direction in degrees, 0 is x+ axis
% X,Y,Z,colp : world data, grasses and color
% hfov : horizontal field of view in degrees, result image range from
%   [-hfov/2,hfov/2]
% resolution : image resolution in degrees/pixel
% Note: close figure when changing world parameters (X,Y,Z,colp), hfov and
%   resolution

z0=0;%getHeight(x,y,X,Y,Z);

%% prepare a figure or reuse figure in previous calls.
% Mike edit. 30/8/13 - always start by sampling the image at 1 degrees 
% resolution and full panorama. This has been hard coded for now. See my 
% later additions for how the code then outputs at the required field of 
% view and resolution.
f = 99;
if ~ishandle(f)
    figure(f);
    pause(0.25);
end
a = get(f,'CurrentAxes');
if isempty(a)
    figure(f);
    a = gca;
    set(f,'Color','c');
    set(f,'Name','Image Grabber');
%     set(f,'Renderer','OpenGL');
    set(f,'Renderer','zbuffer');
    %imgwidth=ceil(hfov/resolution);
    %imgheight = ceil(imgwidth/hfov*75);
    imgwidth=ceil(360);
    imgheight = ceil(imgwidth/360*75);
    set(f, 'Position', [1 500 imgwidth imgheight]);
    set(a,'ActivePositionProperty','position');
    set(a, 'Position', [0 0 1 1]);
end

ph = findobj(get(a,'Children'),'UserData','plothere');
if isempty(ph)
    clear ph
else
    ph2 = findobj(get(a,'Children'),'UserData','plothere2');
    ph3 = findobj(get(a,'Children'),'UserData','plothere3');
end

%% data projection from cartesian coordinates to spherical
[TH,PHI,R]=cart2sph(X-x,Y-y,abs(Z)-z-z0);
% pi to pi: [-180,180)
TH = pi2pi(TH-th/180*pi);

ind=(max(TH')-min(TH')<pi);
A1 = TH(ind,:);
E1 = PHI(ind,:);
D1 = R(ind,:);
c1 = colp(ind,:);

A2 = TH(~ind,:);
E2 = PHI(~ind,:);
D2 = R(~ind,:);
c2 = colp(~ind,:);

A3 = A2;
A3(A3<=0) = A3(A3<=0)+2*pi;
A4 = A2;
A4(A4>0) = A4(A4>0)-2*pi;
%% put data on the figure
if exist('ph','var')
    % change existing data
    set(ph,'XData',D1');
    set(ph,'YData',A1');
    set(ph,'ZData',E1');
    set(ph,'CData',c1');
    
    set(ph2,'XData',D2');
    set(ph2,'YData',A3');
    set(ph2,'ZData',E2');
    set(ph2,'CData',c2');
    
    set(ph3,'XData',D2');
    set(ph3,'YData',A4');
    set(ph3,'ZData',E2');
    set(ph3,'CData',c2');
else
    grasscolormap = zeros(64,3);
    grasscolormap(:,2) = linspace(0,1,64);
    % plot grasses
    ph = patch(D1',A1',E1',c1','EdgeColor','none');
    set(ph,'Userdata','plothere');
    colormap(grasscolormap);
    
    hold on
    ph2 = patch(D2',A3',E2',c2','EdgeColor','none');
    set(ph2,'Userdata','plothere2');
    colormap(grasscolormap);
    
    ph3 = patch(D2',A4',E2',c2','EdgeColor','none');
    set(ph3,'Userdata','plothere3');
    colormap(grasscolormap);
    
    % plot ground
    Xp = [-10 -10 10.5 10.5]'; Yp = [-pi pi pi -pi]';
    Zp = atan2(-z-z0,Xp);
    groundcolor = [229 183 90]/255;
    patch(Xp,Yp,Zp,groundcolor,'EdgeColor','none');
    hold off
    
    % resize and crop
    % note: The view is looking in x+ axis direction.
    %  Thus, horizontal axis is y+ to the left
    %  or (hfov/2) to (-hfov/2) from left to right
    axis equal
    axis off
    %hfov = hfov/180/2*pi;
    %axis([0 14 -hfov hfov -pi/12 pi/3]);
    axis([0 14 -360/180/2*pi 360/180/2*pi -pi/12 pi/3]);
    view([-90 0]);
end

% keyboard;
drawnow
F = getframe(a);
% hgexport(f, 'test.png', options);

img = F.cdata;
% keyboard;
% This is Mike's addition to the code (29/8/13).  The original resize had
% some issues and also the cropping caused. I have deleted all references
% to hfov and resolution from the above.  

% Set image to just the green channel.
img=img(:,:,2);

% Trim the image to the value defined by hfov.
if hfov~=360,
    img=img(:,181-hfov/2:180+hfov/2,:);
end

% Use blockproc command to smooth and subsample the image to that of the
% ant eye
if resolution==1,
    %nothing to be done
elseif resolution==4 && hfov~=360, % Ant eye subsampling
    fun=@(block_struct) uint8(mean(mean(block_struct.data(2:3,2:3))));
    img=blockproc(img,[4 4],fun);
end



function x=pi2pi(x)
x=mod(x,2*pi);
x=x-(x>pi)*2*pi;
% Ai = x>pi;
% x(Ai) = x(Ai)-2*pi;