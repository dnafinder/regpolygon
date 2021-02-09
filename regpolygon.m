function regpolygon(varargin)
%In Euclidean geometry, a regular polygon is a polygon that is equiangular
%(all angles are equal in measure) and equilateral (all sides have the same
%length). This function gives all the possible properties of a regulare
%polygon. Moreover, for polygons with 3<=n<=12 sides, will be shown an
%animation to how contruct it (exactly or approximately) using ruler and
%compass. 
% 
% Syntax: 	regpolygon(varargin)
%    
%     Inputs:
%           n: sides number, a scalar >=3 (3 by default)
%           L: side length, a scalar >0 (1 by default) 
%
%     Outputs:
%           Sides
%           Length
%           Fixed number
%           Apotema
%           Inscribed circle area
%           Height
%           number of diagonals
%           ways of fan triangulations
%           Perimeter
%           Area fixed number (phi)
%           Area
%           Circumradius
%           Circumscribed circle area
%           Interior angle 
%           Exterior angle 
%           Constructible or not 
%
%      Example: 
%
%           Calling on Matlab the function: regpolygon(5)
%
%           Answer is:
%
%                                   Value 
%                                  _______
% 
%     Sides                              5
%     Length                             1
%     Fixed_Number_(f)             0.68819
%     Apotema_(green)              0.68819
%     Inscribed_circle_area         1.4879
%     Height                        1.5388
%     n_of_Diagonals                     5
%     ways_of_fan_triangulations         5
%     Perimeter                          5
%     Area_fixed_number_(phi)       1.7205
%     Area                          1.7205
%     Circumradius_(blue)          0.85065
%     Circumscribed_circle_area     2.2733
% 
% Interior angle: 3/5*pi	108.00°
% Exterior angle: 2/5*pi	72.00°
% Constructible polygon
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo.75@gmail.com
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2020). Regular polygons: facts and drawing
% https://github.com/dnafinder/regpolygon

p = inputParser;
addOptional(p,'n',3,@(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','nonempty','integer','nonzero','>=',3}));
addOptional(p,'L',1,@(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','nonempty','integer','nonzero','>',0}));
parse(p,varargin{:});
n=p.Results.n; L=p.Results.L;
clear p varargin

close all
clc
fn=@(n) 1/2*cot(pi/n);
Apotema=@(n,L) L*fn(n);
Perimeter=@(n,L) n*L;
phi=@(n) n/2*fn(n);
Area=@(n,L) L^2*phi(n);
circle=@(r) pi*r^2;
a=Apotema(n,L);
cr=@(n,L) L/2*csc(pi/n);
circumradius=cr(n,L);
iangle=@(n)(n-2)/n;
eangle=@(n) 2/n;
ndiags=@(n) 1/2*n*(n-3);
if mod(n,2)==0
    h=a*2;
else
    h=a+circumradius;
end
triang=@(n) round(exp(gammaln(2*n+1)-sum(gammaln([n+2 n+1]))));
if n>3
    tri=triang(n-2);
else
    tri=0;
end
    
disp(array2table([n;L;fn(n);a;circle(a);h;ndiags(n);tri;Perimeter(n,L);phi(n);Area(n,L);circumradius;circle(circumradius)],'RowNames',{'Sides','Length','Fixed_Number_(f)','Apotema_(green)','Inscribed_circle_area','Height','n_of_Diagonals','ways_of_fan_triangulations','Perimeter','Area_fixed_number_(phi)','Area','Circumradius_(blue)','Circumscribed_circle_area'},'VariableNames',{'Value'}))
fprintf('Interior angle: %s*pi\t%0.2f°\n',strtrim(rats(iangle(n))),iangle(n)*180);
fprintf('Exterior angle: %s*pi\t%0.2f°\n',strtrim(rats(eangle(n))),eangle(n)*180);
if all(ismember(factor(n),[2 3 5 17 257 65537]))
    disp('Constructible polygon');
else
    disp('Not constructible polygon');
end
clear fn Apotema Perimeter Area circle

if ismember(n,3:1:12)
    gray=repmat(192/255,1,3);
    hold on
    axis equal
    P=zeros(n+1,2); P(2,1)=L;
    %plot the side AB
    plot(P(1:2,1),P(1:2,2),'k','Linewidth',2)
    switch n
        case 3
            %Find perpendicular on mid-point
            [Pa,~]=vesicapiscis([0;0],[L;0]);
            P(3,:)=Pa'; clear Pa
            finishplot(n,L,a,circumradius,P)
        case 4
            P(3:4,:)=[L L; 0 L;];
            %Find perpendicular on mid-point
            vesicapiscis([0;0],[L;0],[8/5*pi,13/5*pi],[7/5*pi,2/5*pi]);
            plot(L/2,0,'kx')
            %Point into mid-point and find the intersection of the arc of
            %radius L/2 and the perpendicular (baricentrum)
            arc(pi,pi/2,L/2,0,L/2,gray,20)
            %Draw the circumscribed circle and find intersection
            ontocircle(n,L,a,circumradius,P)
        case 5
            %Find perpendicular on mid-point
            [Pa,~]=vesicapiscis([0;0],[L;0],[-12/10,2/3*pi],[7/5*pi,pi/3]);
            comet([Pa(1) Pa(1)],[Pa(2) 9/5*L])
            plot([Pa(1) Pa(1)],[Pa(2) 9/5*L],'Color',gray,'Linewidth',2)
            clear Pa
            rextremeperpendicular(L,h,gray)
            %intersection between perpendicular onto right extremity and second
            %arc (BK)
            plot(L,L,'kx')
            comet([L L],[0 L])
            plot([L L],[0 L],'Color',gray,'Linewidth',2)
            %distance MK
            r=L/2*sqrt(5);
            %point into M and draw an arc of radius r
            arc(pi/2,-1/50*pi,L/2,0,r,gray,20)
            %find the intersection between side and ark (AJ)
            [Pa,~]=linecircint([L/2,0],[L,0],[L/2,0],r);
            r=Pa(1); clear Pa
            plot(r,0,'kx')
            comet([L r],[0 0])
            plot([L r],[0 0],'Color',gray,'Linewidth',2)
            %point into A and draw an arc of radius AJ
            arc(0,pi/2,0,0,r,gray,20)
            %find intersection between this arc and the second drawn arc (C)
            [~,Pb]=circcircint([L,0],L,[0,0],r);
            P(3,:)=Pb';
            clear Pb
            plot(P(3,1),P(3,2),'kx')
            %find intersection between this arc and the perpendicular into
            %mid-point (D)
            [g,~]=linecircint([L/2,0],[L/2,L],[0,0],r);
            P(4,:)=[L/2 g(2)];
            plot(P(4,1),P(4,2),'kx')
            clear g
            %distance between B and D
            t=sqrt(L^2/4+P(4,2)^2);
            %point into B and draw an arc of radius t
            arc(pi,pi/2,L,0,t,gray,20)
            %find intersection between this arc and the first drawn arc (E)
            [~,Pb]=circcircint([L,0],t,[0,0],L);
            P(5,:)=Pb'; clear Pb
            plot(P(5,1),P(5,2),'kx')
            finishplot(n,L,a,circumradius,P)
        case 6
            %Find perpendicular on mid-point
            vesicapiscis([0;0],[L;0],[-12/10,2/3*pi],[7/5*pi,pi/3]);
            plot(L/2,a,'ko','Markersize',5,'MarkerFaceColor','k')
            plotcircumscribedcircle(n,L,a,circumradius)
            %find intersections (C)
            [~,Pb]=circcircint([L;0],L,[L/2;a],circumradius);
            P(3,:)=Pb';
            plot(P(3,1),P(3,2),'kx')
            %find intersections (F)
            [Pa,~]=circcircint([0;0],L,[L/2;a],circumradius);
            P(6,:)=Pa';
            plot(P(6,1),P(6,2),'kx')
            clear Pa Pb
            %find intersections (D)
            arc(pi,6*pi/10,P(3,1),P(3,2),L,gray,20)
            [~,Pb]=circcircint(P(3,:),L,[L/2;a],circumradius);
            P(4,:)=Pb'; clear Pb
            plot(P(4,1),P(4,2),'kx')
            %find intersections (E)
            arc(0,2*pi/5,P(6,1),P(6,2),L,gray,20)
            [Pa,~]=circcircint(P(6,:),L,[L/2;a],circumradius);
            P(5,:)=Pa'; clear Pa
            plot(P(5,1),P(5,2),'kx')
            %plot sides
            comet(P(:,1),P(:,2))
            plot(P(:,1),P(:,2),'r','Linewidth',3)
            %plot inscribed circle
            plotinscribedcircle(L,a)
        case 7
            %extend AB by L
            arc(pi,-1/10*pi,L,0,L,gray,20)
            plot(2*L,0,'kx')
            plot([L 2*L],[0 0],'Color',gray,'Linewidth',2)
            %Plot the perpendicular to extremity
            rextremeperpendicular(L,h,gray)
            r=L/4;
            %find the intersection between the perpendicular and the arc
            %obtained pointing in A with radius 2L (K)
            arc(0,pi/2,0,0,2*L,gray,20)
            [Pa,~]=linecircint([L;0],[L;L],[0;0],2*L);
            plot(L,Pa(2),'kx');
            %Plot the angle BAK
            plot([0 L],[0 Pa(2)],'Color',gray,'Linewidth',2)
            %Point into A and draw an arc that you want that intercept BAK
            arc(0,pi/3,0,0,r,gray,10)
            plot(r,0,'kx');
            [Pa,~]=linecircint([0;0],[L;Pa(2)],[0;0],r);
            plot(Pa(1),Pa(2),'kx');
            %Point into intercepts and draw arcs with a radius grater than
            %the half of intercepts distance
            d=(sqrt((Pa(1)-r)^2+Pa(2)^2))*7/8;
            arc(2*pi/3,pi/10,r,0,d,gray,5)
            arc(pi/3,-2/10*pi,Pa(1),Pa(2),d,gray,5)
            [~,Pb]=circcircint([r;0],d,Pa,d);
            clear r d Pa
            plot(Pb(1),Pb(2),'kx');
            r=Pb(2)/Pb(1)*L;
            plot(L,r,'kx');
            plot([0 Pb(1) L],[0 Pb(2) r],'Color',gray,'Linewidth',2)
            d=sqrt(L^2+r^2);
            clear Pb Pa
            %Point into A and B with radius d and you wll have the
            %baricentrum
            arc(pi/6,pi/2,0,0,d,gray,20)
            arc(5/6*pi,pi/2,L,0,d,gray,20)
            ontocircle(n,L,a,circumradius,P)
        case 8
            %Find perpendicular on mid-point
            vesicapiscis([0;0],[L;0]);
            %Plot semicircle in M with radius=AM
            arc(0,pi,L/2,0,L/2,gray,20)
            plot(L/2,L/2,'kx');
            arc(0,2*pi,L/2,L/2,L*sqrt(2)/2,gray,20)
            plot(L/2,L/2*(1+sqrt(2)),'kx');
            ontocircle(n,L,a,circumradius,P)
        case 9
            %Find perpendicular on mid-point
            [Pa,~]=vesicapiscis([0;0],[L;0]);
            plot([Pa(1) Pa(1)],[Pa(2) 11/10*h],'Color',gray,'Linewidth',2)
            %point into intersection and draw and arc of radius L/2 that
            %intercept the perpendicular
            arc(pi/2-4/10,pi/2+4/10,L/2,Pa(2),L/2,gray,10)
            ontocircle(n,L,a,circumradius,P)
        case 10
            %Find perpendicular on mid-point
            vesicapiscis([0;0],[L;0]);
            %Plot the perpendicular to extremity
            rextremeperpendicular(L,h/6,gray)
            %report L/2 onto perpendicular
            arc(pi,pi/2,L,0,L/2,gray,10)
            plot(L,L/2,'kx');
            plot([0 3/2*L],[0 3/4*L],'Color',gray,'Linewidth',2)
            arc(-1/2*pi,2*pi/5,L,L/2,L/2,gray,10)
            [Pa,~]=linecircint([0;0],[3/2*L;3/4*L],[L;L/2],L/2);
            plot(Pa(1),Pa(2),'kx')
            plot([0,Pa(1)],[0 Pa(2)],'Color',gray,'Linewidth',2)
            angle=atan(Pa(2)/Pa(1));
            arc(angle,(angle/pi+2)*pi,0,0,circumradius,[0 0 1],30)
            plot([0 L],[0 0],'Color','r','Linewidth',2)
            plot([circumradius-L circumradius],[0 0],'Color',gray,'Linewidth',2)
            plot(circumradius,0,'kx')
            P(1,:)=[circumradius 0]; 
            angle=2*pi/10;
            for I=2:11
                P(I,:)=circumradius.*[cos(angle*(I-1)) sin(angle*(I-1))];
                plot(P(I,1),P(I,2),'kx')
            end
            %plot sides
            comet(P(:,1),P(:,2))
            plot(P(:,1),P(:,2),'r','Linewidth',3)
            %plot inscribed circle
            plot([0 0],[0 -a],'g','Linewidth',2)
            arc(3/2*pi,7/2*pi,0,0,a,[0 1 0],40)
            plot([0 P(8,1)],[0 P(8,2)],'b','Linewidth',2)
        case 11
            plot([-circumradius circumradius],[0 0],'Color',[0 0 1],'Linewidth',2)
            arc(0,2*pi,0,0,circumradius,[0 0 1],40)
            %%Find perpendicular on mid-point
            vesicapiscis([-circumradius;0],[circumradius;0]);
            plot(0,circumradius,'kx')
            arc(pi,2*pi,0,circumradius,circumradius,gray,20)
            [Pa,Pb]=circcircint([0 0],circumradius,[0 circumradius],circumradius);
            plot(Pa(1),Pa(2),'kx')
            plot(Pb(1),Pb(2),'kx')
            plot([Pa(1) Pb(1)],[Pa(2) Pb(2)],'Color',gray,'Linewidth',2);
            [Pa,~]=vesicapiscis(Pa,[0;Pa(2)]);
            plot(Pa(1),Pb(2),'kx')
            [~,Pb]=linecircint([0;0],[Pa(1);Pb(2)],[0;0],circumradius);
            plot(Pb(1),Pb(2),'ro')
            P(1,:)=Pb'; P(12,:)=P(1,:);
            plot([0 Pb(1)],[0 Pb(2)],'Color',gray,'Linewidth',2);
            r=sqrt((Pb(1)+circumradius)^2+Pb(2)^2);
            arc(pi/2,-pi/2,-circumradius,0,r,gray,20)
            [Pa,~]=circcircint([0 0],circumradius,[-circumradius 0],r);
            P(4,:)=Pa';
            plot(Pa(1),Pa(2),'ro')
            r=2*abs(Pa(2));
            arc(pi/2,0,Pa(1),Pa(2),r,gray,20)
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(7,:)=Pa';
            arc(5*pi/8,9*pi/16,Pa(1),Pa(2),r,gray,10)
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(10,:)=Pa';
            arc(19*pi/16,35*pi/32,Pa(1),Pa(2),r,gray,10)
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(2,:)=Pa';
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(5,:)=Pa';
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(8,:)=Pa';
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(11,:)=Pa';
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(3,:)=Pa';
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(6,:)=Pa';
            [Pa,~]=circcircint([0 0],circumradius,Pa,r);
            plot(Pa(1),Pa(2),'ro')
            P(9,:)=Pa';
            comet(P(:,1),P(:,2))
            plot(P(:,1),P(:,2),'r','Linewidth',3) 
            m=mean(P(8:9,:));
            plot([0 m(1)],[0 m(2)],'g','Linewidth',2)
            m1=P(8,2)/P(8,1); m2=P(9,2)/P(9,1);
            theta=atan(abs((m2-m1)/(1+m1*m2)));
            arc(theta,theta+2*pi,0,0,a,[0 1 0],40)
        case 12
           %Find perpendicular on mid-point
            [Pa,~]=vesicapiscis([0;0],[L;0]);
            plot([Pa(1) Pa(1)],[Pa(2) 3/5*h],'Color',gray,'Linewidth',2)
            r=sqrt(sum(Pa.^2));
            arc(3/2*pi,pi/2,Pa(1),Pa(2),r,gray,20)
            Pb=Pa+[0;r];
            plot(Pb(1),Pb(2),'kx')
            ontocircle(n,L,a,circumradius,P)
    end
    clear gray circumradius a P
    hold off
else
end
end

function arc(start,stop,xc,yc,r,color,points)
xarc=@(r,theta,xc) xc+r.*cos(theta);
yarc=@(r,theta,yc) yc+r.*sin(theta);
theta=linspace(start,stop,points);
comet(xarc(r,theta,xc),yarc(r,theta,yc))
theta=linspace(start,stop,points*3);
plot(xarc(r,theta,xc),yarc(r,theta,yc),'Linewidth',2,'Color',color)
end

function [Pa,Pb]=circcircint(C1,r1,C2,r2)
C1=C1(:); C2=C2(:);
d1=C2-C1;
d2=sum(d1.^2);
P0=(C1+C2)/2+(r1^2-r2^2)/d2/2*d1;
t=((r1+r2)^2-d2)*(d2-(r2-r1)^2);
if t<0
    fprintf('The circles don''t intersect.\n')
else
    T=sqrt(t)/d2/2*[0 -1;1 0]*d1;
    Pa=P0+T;
    Pb=P0-T;
end
end

function [Pa,Pb]=linecircint(P1,P2,C,r)
P1=P1(:); P2=P2(:); C=C(:);
m=(P2(2)-P1(2))/(P2(1)-P1(1)); 
if ~isinf(m)
    q=P1(2)-m*P1(1);
    x=roots([(1+m^2) 2*(m*(q-C(2))-C(1)) C(1)^2+C(2)^2+q^2-2*q*C(2)-r^2]);
    y=m.*x+q;
    Pa=[x(1);y(1)];
    Pb=[x(2);y(2)];
else
    if P1(1)==P2(1) %x=k
        k=P1(1);
        if k^2==r^2
            Pa=[-r;0];
            Pb=[r;0];
        elseif k^2<r^2
            d=sqrt(r^2-k^2);
            Pa=[k;d];
            Pb=[k;-d];
        else
            error('Line can''t intercept circle')
        end
    else
        k=P1(2);
        if k^2==r^2
            Pa=[0;-r];
            Pb=[0;r];
        elseif k^2<r^2
            d=sqrt(r^2-k^2);
            Pa=[d;k];
            Pb=[-d;k];
        else
            error('Line can''t intercept circle')
        end
    end
end
end

function [Pa,Pb]=vesicapiscis(P1,P2,varargin)
p = inputParser;
addRequired(p,'P1');
addRequired(p,'P2');
addOptional(p,'arc1',[8/5*pi,12/5*pi]);
addOptional(p,'arc2',[7/5*pi,3/5*pi]);
addOptional(p,'color',repmat(192/255,1,3));
parse(p,P1,P2,varargin{:});
arc1=p.Results.arc1; arc2=p.Results.arc2;
color=p.Results.color; 
clear p varargin

%Find perpendicular on mid-point
L=sqrt((P2(1)-P1(1))^2+(P2(2)-P1(2))^2);
arc(arc1(1),arc1(2),P1(1),P1(2),L,color,20)
arc(arc2(1),arc2(2),P2(1),P2(2),L,color,20)
%find intersection
[Pa,Pb]=circcircint(P1,L,P2,L);
plot(Pa(1),Pa(2),'kx')
plot(Pb(1),Pb(2),'kx')
plot([Pb(1) Pa(1)],[Pb(2) Pa(2)],'Color',color,'Linewidth',2)
end

function rextremeperpendicular(L,h,color)
%Plot the perpendicular to the right extremity
r=L/4; %choose a radius
%Plot a semicircle pointing into B with radius r
arc(pi,-pi/10,L,0,r,color,20)
%extend the side until it intercept the semicircle
plot(L+r,0,'kx')
plot([L L+r],[0 0],'Color',color,'Linewidth',2)
%Find the mid-point
[Pa,~]=vesicapiscis([L-r;0],[L+r;0]);
plot([L L],[Pa(2) 11/10*h],'Color',color,'Linewidth',2)
end

function finishplot(n,L,a,circumradius,P)
%plot sides
comet(P(:,1),P(:,2))
plot(P(:,1),P(:,2),'r','Linewidth',3)
plot(L/2,a,'ko','Markersize',5,'MarkerFaceColor','k')
%plot inscribed circle
plotinscribedcircle(L,a)
%plot circumscribed circle
plotcircumscribedcircle(n,L,a,circumradius)
end

function plotinscribedcircle(L,a)
%plot inscribed circle
plot([L/2 L/2],[a 0],'g','Linewidth',2)
arc(3/2*pi,7/2*pi,L/2,a,a,[0 1 0],40)
end

function plotcircumscribedcircle(n,L,a,circumradius)
%plot circumscribed circle
plot([0 L/2],[0 a],'b','Linewidth',2)
angle=(pi*(n-2))/n/2+pi;
arc(angle,angle+2*pi,L/2,a,circumradius,[0 0 1],40)
end

function ontocircle(n,L,a,circumradius,P)
plot(L/2,a,'ko','Markersize',5,'MarkerFaceColor','k')
plotcircumscribedcircle(n,L,a,circumradius)
%find intersections
for I=3:n
    [~,Pb]=circcircint(P(I-1,:),L,[L/2;a],circumradius);
    P(I,:)=Pb';
    plot(P(I,1),P(I,2),'kx')
end
%plot sides
comet(P(:,1),P(:,2))
plot(P(:,1),P(:,2),'r','Linewidth',3)
%plot inscribed circle
plotinscribedcircle(L,a)
end
