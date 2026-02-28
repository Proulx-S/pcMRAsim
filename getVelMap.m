function v = getVelMap(r, ID, type, plugDiam)
% r        : radial coordinate from the center of the grid [mm]
% ID       : Inner diameter of the vessel [mm]
% type     : [type num2str(p2op1)] type of velocity profile ('parabolic')
% p2op1    : [0 to 1] 1 is pure quadratic, 0 is pure linear, inbetween is a mixture
% plugDiam : Diameter of the central plug flow region [mm]
% v        : Velocity, arbitrarily scaled [cm/s]

if contains(type, 'parabolic')
    p2op1 = replace(type, 'parabolic','');
    if isempty(p2op1)
        p2op1 = 1; 
    else
        p2op1 = str2double(p2op1);
    end
    type = 'parabolic';
end
if ~exist('plugDiam', 'var'); plugDiam = []; end
if isempty(plugDiam        ); plugDiam = 0 ; end

plugR = plugDiam/2;
lamiR = (ID-plugDiam)/2; %
rLami = r./lamiR - plugR;


switch type
    case 'parabolic'

        % linear/parabolic flow profile
        p0 = 1              ;
        p1 = 1 .* -(1-p2op1);
        p2 = 1 .*    -p2op1 ;
        v = p2.*((r-plugR)./(lamiR)).^2 + p1.*((r-plugR)./(lamiR)) + p0; % full quadratic flow profile
        v(v<0|isnan(v)) = 0;
        % plug flow center
        v(r<=plugR) = 1;

    % case 'plugFlow'
    %     v = zeros(size(r));
    %     v(r<ID/2) = 1; % 1 inside, 0 outside
    % case 'bluntedParabolic'
    %     v = 1 - (r/(ID/2)).^2; % parabola (vMax at the center, 0 at the edge, negative beyond)
    %     v(v<0) = 0; % make it a cross-section profile (0 beyond the edge)
    %     v(r<plugDiam/2) = 1 - ((plugDiam/2)/(ID/2)).^2; % set the plug flow region to the same velocity as the parabolic profile at the plug diameter
    otherwise
        error('Invalid velocity profile type');
end