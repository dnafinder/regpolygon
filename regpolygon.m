function S = regpolygon(varargin)
%REGPOLYGON Properties and construction of a regular polygon.
%
%   Syntax
%   ------
%   regpolygon
%   regpolygon(n)
%   regpolygon(n, L)
%   S = regpolygon(...)
%
%   regpolygon(n, L, 'Animate', false, 'CloseFigures', false, 'ClearCommandWindow', false)
%
%   Description
%   -----------
%   REGPOLYGON computes several geometric properties of a regular polygon
%   in Euclidean geometry. A regular polygon is equiangular (all interior
%   angles equal) and equilateral (all sides have the same length).
%
%   The function:
%     - computes side-related, area-related, and circle-related quantities;
%     - reports interior/exterior angles and constructibility;
%     - displays the results in a MATLAB table in the Command Window;
%     - for polygons with 3 <= n <= 12, shows an animation of how to
%       construct the polygon (exactly or approximately) using ruler and
%       compass.
%
%   Inputs
%   ------
%   Positional inputs:
%     regpolygon()
%         Uses n = 3, L = 1.
%
%     regpolygon(n)
%         Uses side length L = 1.
%
%     regpolygon(n, L)
%         Uses user-defined number of sides and side length.
%
%     n : Number of sides of the regular polygon.
%         Type: scalar numeric, integer, real, finite, >= 3.
%         Default: 3.
%
%     L : Side length of the polygon.
%         Type: scalar numeric, real, finite, > 0.
%         Default: 1.
%
%   Name-Value options:
%     'ClearCommandWindow' : Logical or 0/1 flag that controls whether the
%                            Command Window is cleared at the beginning.
%                            Default: true.
%
%     'CloseFigures'       : Logical or 0/1 flag that controls whether all
%                            open figures are closed at the beginning.
%                            Default: true.
%
%     'Animate'            : Logical or 0/1 flag that controls whether the
%                            geometric construction is animated (comet
%                            traces and incremental drawing).
%                            true  (1) - full animation (original behavior).
%                            false (0) - draw only the final lines and
%                                        circles without comet-style
%                                        motion.
%                            Default: true.
%
%   Outputs
%   -------
%   By default (when called without output), the function:
%     - prints a table in the Command Window with the following rows:
%
%       Sides
%       Length
%       Fixed_Number_(f)
%       Apotema_(green)
%       Inscribed_circle_area
%       Height
%       n_of_Diagonals
%       ways_of_fan_triangulations
%       Perimeter
%       Area_fixed_number_(phi)
%       Area
%       Circumradius_(blue)
%       Circumscribed_circle_area
%
%     - prints interior and exterior angles in both radians (as rational
%       multiples of pi) and degrees;
%     - prints whether the polygon is constructible or not.
%
%   When an output is requested:
%     S = regpolygon(...)
%     returns a structure S with fields:
%       S.Sides
%       S.Length
%       S.FixedNumber_f
%       S.Apothem
%       S.InscribedCircleArea
%       S.Height
%       S.NumDiagonals
%       S.FanTriangulations
%       S.Perimeter
%       S.AreaPhi
%       S.Area
%       S.Circumradius
%       S.CircumscribedCircleArea
%       S.InteriorAngleFactorPi   (fraction of pi)
%       S.ExteriorAngleFactorPi   (fraction of pi)
%       S.InteriorAngleDeg
%       S.ExteriorAngleDeg
%       S.Constructible           (logical)
%
%   Example
%   -------
%   % Regular pentagon with unit side length
%   regpolygon(5)
%
%   % Example output (abbreviated):
%   %
%   %                                 Value
%   %                                _______
%   %   Sides                             5
%   %   Length                            1
%   %   Fixed_Number_(f)             0.68819
%   %   Apotema_(green)              0.68819
%   %   Inscribed_circle_area         1.4879
%   %   ...
%   %
%   %   Interior angle: 3/5*pi    108.00°
%   %   Exterior angle: 2/5*pi     72.00°
%   %   Constructible polygon
%
%   % Obtain the properties programmatically:
%   S = regpolygon(5, 1);
%
%   Notes
%   -----
%   - For 3 <= n <= 12, an animated geometric construction of the regular
%     polygon is displayed, using ruler-and-compass-style steps. Set
%     'Animate' to false to skip the comet-based animation and draw only
%     the final figure.
%   - For n outside this range, only the numerical properties are shown.
%   - The constructibility test checks whether n has only Fermat primes
%     and powers of 2 as prime factors (Gauss–Wantzel theorem). For very
%     large n, the factorization step may be time-consuming.
%   - 'CloseFigures' and 'ClearCommandWindow' allow you to embed this
%     function into larger scripts without resetting your graphical or
%     console state.
%
%   Citation
%   --------
%   If you use this function in academic or technical work, please cite:
%
%     Cardillo G. (2020)
%     "Regular polygons: facts and drawing".
%     Available from GitHub:
%     https://github.com/dnafinder/regpolygon
%
%   Metadata
%   --------
%   Author : Giuseppe Cardillo
%   Email  : giuseppe.cardillo.75@gmail.com
%   GitHub : https://github.com/dnafinder
%   Created: 2020-01-01
%   Updated: 2025-11-20
%   Version: 2.1.0
%
%   License
%   -------
%   This function is distributed under the MIT License.
%   See the LICENSE file in the GitHub repository for details.
%

% ---------------------------
% Input parsing and checking
% ---------------------------
p = inputParser;

% Positional syntax: regpolygon(n, L)
addOptional(p,'n',3,@(x) validateattributes(x,{'numeric'}, ...
    {'scalar','real','finite','nonnan','nonempty','integer','>=',3}));
addOptional(p,'L',1,@(x) validateattributes(x,{'numeric'}, ...
    {'scalar','real','finite','nonnan','nonempty','>',0}));

% Behavior controls
addParameter(p,'ClearCommandWindow',true,@(x) islogical(x) || ...
    (isnumeric(x) && isscalar(x) && any(x==[0 1])));
addParameter(p,'CloseFigures',true,@(x) islogical(x) || ...
    (isnumeric(x) && isscalar(x) && any(x==[0 1])));
addParameter(p,'Animate',true,@(x) islogical(x) || ...
    (isnumeric(x) && isscalar(x) && any(x==[0 1])));

parse(p,varargin{:});
n         = p.Results.n;
L         = p.Results.L;
doClc     = logical(p.Results.ClearCommandWindow);
doClose   = logical(p.Results.CloseFigures);
doAnimate = logical(p.Results.Animate);

clear p varargin

% Store animation flag for helper functions (file-local)
regpolygonAnimateFlag(doAnimate);

% ---------------------------
% Console and figure handling
% ---------------------------
if doClose
    close all
end
if doClc
    clc
end

% ---------------------------
% Geometric helper functions
% ---------------------------
% Fixed number f related to area and apothem
fn        = @(n) 0.5 * cot(pi/n);
Apotema   = @(n,L) L * fn(n);
Perimeter = @(n,L) n * L;
phi       = @(n) n/2 * fn(n);         % area fixed number phi
Area      = @(n,L) L^2 * phi(n);
circle    = @(r) pi * r^2;            % area of a circle of radius r

% Apothem and circumradius
a            = Apotema(n,L);
cr           = @(n,L) L/2 * csc(pi/n);
circumradius = cr(n,L);

% Angles (as fractions of pi)
iangle = @(n) (n-2)/n;   % interior angle / pi
eangle = @(n) 2/n;       % exterior angle / pi

% Number of diagonals
ndiags = @(n) 0.5 * n * (n-3);

% Height of the polygon (distance between two parallel sides through center)
if mod(n,2) == 0
    h = 2 * a;
else
    h = a + circumradius;
end

% Ways of fan triangulations (from one fixed vertex) via Catalan number C_k
% For an n-gon, the number is Catalan(n-2).
triang = @(k) nchoosek(2*k, k) / (k+1);
if n > 3
    tri = triang(n-2);
else
    tri = 0;
end

% ---------------------------
% Build and display result table
% ---------------------------
valSides        = n;
valLength       = L;
valFixed        = fn(n);
valApothem      = a;
valInscribed    = circle(a);
valHeight       = h;
valDiagonals    = ndiags(n);
valTriang       = tri;
valPerimeter    = Perimeter(n,L);
valPhi          = phi(n);
valArea         = Area(n,L);
valCircumradius = circumradius;
valCircArea     = circle(circumradius);

T = array2table( ...
    [ valSides; ...
      valLength; ...
      valFixed; ...
      valApothem; ...
      valInscribed; ...
      valHeight; ...
      valDiagonals; ...
      valTriang; ...
      valPerimeter; ...
      valPhi; ...
      valArea; ...
      valCircumradius; ...
      valCircArea ], ...
    'RowNames', { ...
      'Sides', ...
      'Length', ...
      'Fixed_Number_(f)', ...
      'Apotema_(green)', ...
      'Inscribed_circle_area', ...
      'Height', ...
      'n_of_Diagonals', ...
      'ways_of_fan_triangulations', ...
      'Perimeter', ...
      'Area_fixed_number_(phi)', ...
      'Area', ...
      'Circumradius_(blue)', ...
      'Circumscribed_circle_area'}, ...
    'VariableNames', {'Value'});

disp(T);

% Interior and exterior angles in radians (multiple of pi) and degrees
intFactor = iangle(n);
extFactor = eangle(n);
intDeg    = intFactor * 180;
extDeg    = extFactor * 180;

fprintf('Interior angle: %s*pi\t%0.2f°\n', ...
    strtrim(rats(intFactor)), intDeg);
fprintf('Exterior angle: %s*pi\t%0.2f°\n', ...
    strtrim(rats(extFactor)), extDeg);

% Constructibility test (Gauss–Wantzel)
factors        = factor(n);
isConstructible = all(ismember(factors, [2 3 5 17 257 65537]));
if isConstructible
    disp('Constructible polygon');
else
    disp('Not constructible polygon');
end

% ---------------------------
% Ruler-and-compass animation
% ---------------------------
if ismember(n,3:1:12)
    gray = repmat(192/255,1,3);
    hold on
    axis equal
    animate = regpolygonAnimateFlag(); % current animation mode

    % P will collect the vertices of the polygon as 2D points
    P       = zeros(n+1,2);
    P(2,1)  = L;

    % Plot the initial side AB on the x-axis
    plot(P(1:2,1), P(1:2,2), 'k', 'Linewidth', 2)

    switch n
        case 3
            % Equilateral triangle construction
            [Pa,~] = vesicapiscis([0;0],[L;0]);
            P(3,:) = Pa'; clear Pa
            finishplot(n,L,a,circumradius,P)

        case 4
            % Square construction
            P(3:4,:) = [L L; 0 L];
            vesicapiscis([0;0],[L;0],[8/5*pi,13/5*pi],[7/5*pi,2/5*pi]);
            plot(L/2,0,'kx')
            arc(pi,pi/2,L/2,0,L/2,gray,20)
            ontocircle(n,L,a,circumradius,P)

        case 5
            % Regular pentagon construction
            [Pa,~] = vesicapiscis([0;0],[L;0],[-12/10,2/3*pi],[7/5*pi,pi/3]);
            if animate
                comet([Pa(1) Pa(1)],[Pa(2) 9/5*L])
            end
            plot([Pa(1) Pa(1)],[Pa(2) 9/5*L],'Color',gray,'Linewidth',2)
            clear Pa
            rextremeperpendicular(L,h,gray)
            plot(L,L,'kx')
            if animate
                comet([L L],[0 L])
            end
            plot([L L],[0 L],'Color',gray,'Linewidth',2)
            r = L/2 * sqrt(5);
            arc(pi/2,-1/50*pi,L/2,0,r,gray,20)
            [Pa,~] = linecircint([L/2,0],[L,0],[L/2,0],r);
            r = Pa(1); clear Pa
            plot(r,0,'kx')
            if animate
                comet([L r],[0 0])
            end
            plot([L r],[0 0],'Color',gray,'Linewidth',2)
            arc(0,pi/2,0,0,r,gray,20)
            [~,Pb] = circcircint([L,0],L,[0,0],r);
            P(3,:) = Pb'; clear Pb
            plot(P(3,1),P(3,2),'kx')
            [g,~] = linecircint([L/2,0],[L/2,L],[0,0],r);
            P(4,:) = [L/2 g(2)];
            plot(P(4,1),P(4,2),'kx')
            clear g
            t = sqrt(L^2/4 + P(4,2)^2);
            arc(pi,pi/2,L,0,t,gray,20)
            [~,Pb] = circcircint([L,0],t,[0,0],L);
            P(5,:) = Pb'; clear Pb
            plot(P(5,1),P(5,2),'kx')
            finishplot(n,L,a,circumradius,P)

        case 6
            % Regular hexagon construction
            vesicapiscis([0;0],[L;0],[-12/10,2/3*pi],[7/5*pi,pi/3]);
            plot(L/2,a,'ko','Markersize',5,'MarkerFaceColor','k')
            plotcircumscribedcircle(n,L,a,circumradius)
            [~,Pb] = circcircint([L;0],L,[L/2;a],circumradius);
            P(3,:) = Pb';
            plot(P(3,1),P(3,2),'kx')
            [Pa,~] = circcircint([0;0],L,[L/2;a],circumradius);
            P(6,:) = Pa';
            plot(P(6,1),P(6,2),'kx')
            clear Pa Pb
            arc(pi,6*pi/10,P(3,1),P(3,2),L,gray,20)
            [~,Pb] = circcircint(P(3,:),L,[L/2;a],circumradius);
            P(4,:) = Pb'; clear Pb
            plot(P(4,1),P(4,2),'kx')
            arc(0,2*pi/5,P(6,1),P(6,2),L,gray,20)
            [Pa,~] = circcircint(P(6,:),L,[L/2;a],circumradius);
            P(5,:) = Pa'; clear Pa
            plot(P(5,1),P(5,2),'kx')
            if animate
                comet(P(:,1),P(:,2))
            end
            plot(P(:,1),P(:,2),'r','Linewidth',3)
            plotinscribedcircle(L,a)

        case 7
            % Regular heptagon construction (approximate)
            arc(pi,-1/10*pi,L,0,L,gray,20)
            plot(2*L,0,'kx')
            plot([L 2*L],[0 0],'Color',gray,'Linewidth',2)
            rextremeperpendicular(L,h,gray)
            r = L/4;
            arc(0,pi/2,0,0,2*L,gray,20)
            [Pa,~] = linecircint([L;0],[L;L],[0;0],2*L);
            plot(L,Pa(2),'kx');
            plot([0 L],[0 Pa(2)],'Color',gray,'Linewidth',2)
            arc(0,pi/3,0,0,r,gray,10)
            plot(r,0,'kx');
            [Pa,~] = linecircint([0;0],[L;Pa(2)],[0;0],r);
            plot(Pa(1),Pa(2),'kx');
            d = (sqrt((Pa(1)-r)^2 + Pa(2)^2)) * 7/8;
            arc(2*pi/3,pi/10,r,0,d,gray,5)
            arc(pi/3,-2/10*pi,Pa(1),Pa(2),d,gray,5)
            [~,Pb] = circcircint([r;0],d,Pa,d);
            clear r d Pa
            plot(Pb(1),Pb(2),'kx');
            r = Pb(2)/Pb(1) * L;
            plot(L,r,'kx');
            plot([0 Pb(1) L],[0 Pb(2) r],'Color',gray,'Linewidth',2)
            d = sqrt(L^2 + r^2);
            clear Pb Pa
            arc(pi/6,pi/2,0,0,d,gray,20)
            arc(5/6*pi,pi/2,L,0,d,gray,20)
            ontocircle(n,L,a,circumradius,P)

        case 8
            % Regular octagon construction
            vesicapiscis([0;0],[L;0]);
            arc(0,pi,L/2,0,L/2,gray,20)
            plot(L/2,L/2,'kx');
            arc(0,2*pi,L/2,L/2,L*sqrt(2)/2,gray,20)
            plot(L/2,L/2*(1+sqrt(2)),'kx');
            ontocircle(n,L,a,circumradius,P)

        case 9
            % Regular nonagon construction (approximate)
            [Pa,~] = vesicapiscis([0;0],[L;0]);
            plot([Pa(1) Pa(1)],[Pa(2) 11/10*h],'Color',gray,'Linewidth',2)
            arc(pi/2-4/10,pi/2+4/10,L/2,Pa(2),L/2,gray,10)
            ontocircle(n,L,a,circumradius,P)

        case 10
            % Regular decagon construction
            vesicapiscis([0;0],[L;0]);
            rextremeperpendicular(L,h/6,gray)
            arc(pi,pi/2,L,0,L/2,gray,10)
            plot(L,L/2,'kx');
            plot([0 3/2*L],[0 3/4*L],'Color',gray,'Linewidth',2)
            arc(-1/2*pi,2*pi/5,L,L/2,L/2,gray,10)
            [Pa,~] = linecircint([0;0],[3/2*L;3/4*L],[L;L/2],L/2);
            plot(Pa(1),Pa(2),'kx')
            plot([0,Pa(1)],[0 Pa(2)],'Color',gray,'Linewidth',2)
            angle = atan(Pa(2)/Pa(1));
            arc(angle,(angle/pi+2)*pi,0,0,circumradius,[0 0 1],30)
            plot([0 L],[0 0],'Color','r','Linewidth',2)
            plot([circumradius-L circumradius],[0 0],'Color',gray,'Linewidth',2)
            plot(circumradius,0,'kx')
            P(1,:) = [circumradius 0];
            angle = 2*pi/10;
            for I = 2:11
                P(I,:) = circumradius .* [cos(angle*(I-1)) sin(angle*(I-1))];
                plot(P(I,1),P(I,2),'kx')
            end
            if animate
                comet(P(:,1),P(:,2))
            end
            plot(P(:,1),P(:,2),'r','Linewidth',3)
            plot([0 0],[0 -a],'g','Linewidth',2)
            arc(3/2*pi,7/2*pi,0,0,a,[0 1 0],40)
            plot([0 P(8,1)],[0 P(8,2)],'b','Linewidth',2)

        case 11
            % Regular hendecagon construction (approximate)
            plot([-circumradius circumradius],[0 0],'Color',[0 0 1],'Linewidth',2)
            arc(0,2*pi,0,0,circumradius,[0 0 1],40)
            vesicapiscis([-circumradius;0],[circumradius;0]);
            plot(0,circumradius,'kx')
            arc(pi,2*pi,0,circumradius,circumradius,gray,20)
            [Pa,Pb] = circcircint([0 0],circumradius,[0 circumradius],circumradius);
            plot(Pa(1),Pa(2),'kx')
            plot(Pb(1),Pb(2),'kx')
            plot([Pa(1) Pb(1)],[Pa(2) Pb(2)],'Color',gray,'Linewidth',2);
            [Pa,~] = vesicapiscis(Pa,[0;Pa(2)]);
            plot(Pa(1),Pb(2),'kx')
            [~,Pb] = linecircint([0;0],[Pa(1);Pb(2)],[0;0],circumradius);
            plot(Pb(1),Pb(2),'ro')
            P(1,:) = Pb'; P(12,:) = P(1,:);
            plot([0 Pb(1)],[0 Pb(2)],'Color',gray,'Linewidth',2);
            r = sqrt((Pb(1)+circumradius)^2 + Pb(2)^2);
            arc(pi/2,-pi/2,-circumradius,0,r,gray,20)
            [Pa,~] = circcircint([0 0],circumradius,[-circumradius 0],r);
            P(4,:) = Pa';
            plot(Pa(1),Pa(2),'ro')
            r = 2 * abs(Pa(2));
            arc(pi/2,0,Pa(1),Pa(2),r,gray,20)
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(7,:) = Pa';
            arc(5*pi/8,9*pi/16,Pa(1),Pa(2),r,gray,10)
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(10,:) = Pa';
            arc(19*pi/16,35*pi/32,Pa(1),Pa(2),r,gray,10)
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(2,:) = Pa';
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(5,:) = Pa';
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(8,:) = Pa';
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(11,:) = Pa';
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(3,:) = Pa';
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(6,:) = Pa';
            [Pa,~] = circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(9,:) = Pa';
            if animate
                comet(P(:,1),P(:,2))
            end
            plot(P(:,1),P(:,2),'r','Linewidth',3)
            m = mean(P(8:9,:));
            plot([0 m(1)],[0 m(2)],'g','Linewidth',2)
            m1 = P(8,2)/P(8,1); m2 = P(9,2)/P(9,1);
            theta = atan(abs((m2-m1)/(1+m1*m2)));
            arc(theta,theta+2*pi,0,0,a,[0 1 0],40)

        case 12
            % Regular dodecagon construction
            [Pa,~] = vesicapiscis([0;0],[L;0]);
            plot([Pa(1) Pa(1)],[Pa(2) 3/5*h],'Color',gray,'Linewidth',2)
            r = sqrt(sum(Pa.^2));
            arc(3/2*pi,pi/2,Pa(1),Pa(2),r,gray,20)
            Pb = Pa + [0;r];
            plot(Pb(1),Pb(2),'kx')
            ontocircle(n,L,a,circumradius,P)
    end

    clear gray circumradius a P
    hold off
end

% ---------------------------
% Optional structured output
% ---------------------------
if nargout > 0
    S = struct();
    S.Sides                 = valSides;
    S.Length                = valLength;
    S.FixedNumber_f         = valFixed;
    S.Apothem               = valApothem;
    S.InscribedCircleArea   = valInscribed;
    S.Height                = valHeight;
    S.NumDiagonals          = valDiagonals;
    S.FanTriangulations     = valTriang;
    S.Perimeter             = valPerimeter;
    S.AreaPhi               = valPhi;
    S.Area                  = valArea;
    S.Circumradius          = valCircumradius;
    S.CircumscribedCircleArea = valCircArea;
    S.InteriorAngleFactorPi = intFactor;
    S.ExteriorAngleFactorPi = extFactor;
    S.InteriorAngleDeg      = intDeg;
    S.ExteriorAngleDeg      = extDeg;
    S.Constructible         = isConstructible;
end

end % regpolygon main function

% -------------------------------------------------------------------------
% Helper: global animation flag (file-local persistent)
% -------------------------------------------------------------------------
function tf = regpolygonAnimateFlag(tfIn)
%REGPOLYGONANIMATEFLAG Get or set animation flag for helper functions.
persistent flag
if nargin > 0
    flag = logical(tfIn);
end
if isempty(flag)
    flag = true;
end
tf = flag;
end

% -------------------------------------------------------------------------
% Helper functions (geometry and drawing)
% -------------------------------------------------------------------------
function arc(start,stop,xc,yc,r,color,points)
%ARC Draw an arc with optional comet-style animation.
animate = regpolygonAnimateFlag();
xarc   = @(r,theta,xc) xc + r .* cos(theta);
yarc   = @(r,theta,yc) yc + r .* sin(theta);
theta  = linspace(start,stop,points);
if animate
    comet(xarc(r,theta,xc),yarc(r,theta,yc))
end
theta  = linspace(start,stop,points*3);
plot(xarc(r,theta,xc),yarc(r,theta,yc),'Linewidth',2,'Color',color)
end

function [Pa,Pb] = circcircint(C1,r1,C2,r2)
%CIRCCIRCINT Intersection between two circles (if any).
C1 = C1(:); C2 = C2(:);
d1 = C2 - C1;
d2 = sum(d1.^2);
P0 = (C1 + C2)/2 + (r1^2 - r2^2)/d2/2 * d1;
t  = ((r1 + r2)^2 - d2) * (d2 - (r2 - r1)^2);
if t < 0
    fprintf('The circles don''t intersect.\n')
else
    T  = sqrt(t)/d2/2 * [0 -1; 1 0] * d1;
    Pa = P0 + T;
    Pb = P0 - T;
end
end

function [Pa,Pb] = linecircint(P1,P2,C,r)
%LINECIRCINT Intersection between a line (through P1,P2) and a circle.
P1 = P1(:); P2 = P2(:); C = C(:);
m  = (P2(2)-P1(2))/(P2(1)-P1(1)); 
if ~isinf(m)
    % Non-vertical line
    q = P1(2) - m * P1(1);
    x = roots([(1+m^2) 2*(m*(q-C(2))-C(1)) C(1)^2 + C(2)^2 + q^2 - 2*q*C(2) - r^2]);
    y = m .* x + q;
    Pa = [x(1); y(1)];
    Pb = [x(2); y(2)];
else
    % Vertical or horizontal special cases
    if P1(1) == P2(1) % x = k
        k = P1(1);
        if k^2 == r^2
            Pa = [-r; 0];
            Pb = [ r; 0];
        elseif k^2 < r^2
            d  = sqrt(r^2 - k^2);
            Pa = [k; d];
            Pb = [k; -d];
        else
            error('Line can''t intercept circle')
        end
    else
        k = P1(2);
        if k^2 == r^2
            Pa = [0; -r];
            Pb = [0;  r];
        elseif k^2 < r^2
            d  = sqrt(r^2 - k^2);
            Pa = [d; k];
            Pb = [-d; k];
        else
            error('Line can''t intercept circle')
        end
    end
end
end

function [Pa,Pb] = vesicapiscis(P1,P2,varargin)
%VESICAPISCIS Classical vesica piscis construction between two points.
p = inputParser;
addRequired(p,'P1');
addRequired(p,'P2');
addOptional(p,'arc1',[8/5*pi,12/5*pi]);
addOptional(p,'arc2',[7/5*pi,3/5*pi]);
addOptional(p,'color',repmat(192/255,1,3));
parse(p,P1,P2,varargin{:});
arc1  = p.Results.arc1;
arc2  = p.Results.arc2;
color = p.Results.color;
clear p varargin

L = sqrt((P2(1)-P1(1))^2 + (P2(2)-P1(2))^2);
arc(arc1(1),arc1(2),P1(1),P1(2),L,color,20)
arc(arc2(1),arc2(2),P2(1),P2(2),L,color,20)
[Pa,Pb] = circcircint(P1,L,P2,L);
plot(Pa(1),Pa(2),'kx')
plot(Pb(1),Pb(2),'kx')
plot([Pb(1) Pa(1)],[Pb(2) Pa(2)],'Color',color,'Linewidth',2)
end

function rextremeperpendicular(L,h,color)
%REXTREMEPERPENDICULAR Perpendicular from the right extremity of a segment.
r = L/4; % choose a radius
arc(pi,-pi/10,L,0,r,color,20)
plot(L+r,0,'kx')
plot([L L+r],[0 0],'Color',color,'Linewidth',2)
[Pa,~] = vesicapiscis([L-r;0],[L+r;0]);
plot([L L],[Pa(2) 11/10*h],'Color',color,'Linewidth',2)
end

function finishplot(n,L,a,circumradius,P)
%FINISHPLOT Final drawing for basic constructions (3,4,5).
animate = regpolygonAnimateFlag();
if animate
    comet(P(:,1),P(:,2))
end
plot(P(:,1),P(:,2),'r','Linewidth',3)
plot(L/2,a,'ko','Markersize',5,'MarkerFaceColor','k')
plotinscribedcircle(L,a)
plotcircumscribedcircle(n,L,a,circumradius)
end

function plotinscribedcircle(L,a)
%PLOTINSCRIBEDCIRCLE Draw the inscribed circle and its radius.
plot([L/2 L/2],[a 0],'g','Linewidth',2)
arc(3/2*pi,7/2*pi,L/2,a,a,[0 1 0],40)
end

function plotcircumscribedcircle(n,L,a,circumradius)
%PLOTCIRCUMSCRIBEDCIRCLE Draw the circumscribed circle and a radius.
plot([0 L/2],[0 a],'b','Linewidth',2)
angle = (pi*(n-2))/n/2 + pi;
arc(angle,angle+2*pi,L/2,a,circumradius,[0 0 1],40)
end

function ontocircle(n,L,a,circumradius,P)
%ONTOCIRCLE Place vertices on circumscribed circle and finish drawing.
animate = regpolygonAnimateFlag();
plot(L/2,a,'ko','Markersize',5,'MarkerFaceColor','k')
plotcircumscribedcircle(n,L,a,circumradius)
for I = 3:n
    [~,Pb] = circcircint(P(I-1,:),L,[L/2;a],circumradius);
    P(I,:) = Pb';
    plot(P(I,1),P(I,2),'kx')
end
if animate
    comet(P(:,1),P(:,2))
end
plot(P(:,1),P(:,2),'r','Linewidth',3)
plotinscribedcircle(L,a)
end
