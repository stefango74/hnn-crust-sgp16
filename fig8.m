% Copyright 2016 Stefan Ohrhallinger
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function main()

% variables are:
% ratio ||x,s1||/||y,s1|| for 1-reach in ]0, 1]
% onaxis = 0 (roll on lfs(y)), onaxis = 1 (roll on x-axis = curve approximation for angle_x < 1)
% angle_[x|y] in [0, 1] where 1 is max. angle with tangent at s1 = horizontal axis (~27 degrees for rho=0.9)

reach(1, 0, 1, 1, 'figure8a');
reach(1/3, 0, 1, 1, 'figure8b');
reach(1/sqrt(2), 0, 1, 1, 'figure8c'); % worst case ratio
reach(1/sqrt(2), 0, 1, 0, 'figure8d');
reach(1/sqrt(2), 1, 0.5, 0, 'figure8e');
reach(1/sqrt(2), 1, 0, 0, 'figure8f');
reach(0.001, 1, 0, 0, 'figure8g'); % avoid exact 0
reach(1, 1, 0, 1, 'figure8h');

end

% show limits for reach sampling
function reach(ratio, onaxis, angle_x, angle_y, name)

% constants
RHO=0.9;
LFSY = 1;
CIRCLECOUNT = 10;
STEPSIZE = 5;

initializeFigure();

[c, s1, yy, s2, s1s2] = computeConstantItems(LFSY, RHO, angle_y);
[xx, z, startH, endH, posH] = computeVariableItems(LFSY, RHO, ratio, angle_x, c, s1, s1s2);

drawLfsDisks(LFSY, xx, c, ratio, CIRCLECOUNT, 'Rx', 0);
drawRtDisks(LFSY, LFSY, c, xx, s1, 'Rts1', STEPSIZE, onaxis, 0);
drawLfsDisks(LFSY, xx, c, ratio, CIRCLECOUNT, 'Rx', 1);
drawRtDisks(LFSY, LFSY, c, xx, s1, 'Rts1', STEPSIZE, onaxis, 1);
disktrans(s1(1), s1(2), s1s2, [1 0 0]);

circle(s1(1), s1(2), LFSY, 'black');
drawItems(LFSY, ratio, c, xx, s1, yy, s2, z, s1s2, startH, endH, posH);

% lfs(y)/2 distance
line([s1(1) yy(1)], [s1(2) yy(2)], 'Color', 'black');

saveFigure(name);

end

function initializeFigure()

fig = figure();
hold on;
fig.ToolBar = 'none';
fig.MenuBar = 'none';
fig.Color = 'white';
fig.Position = [0 0 600 600];
axis manual;
axis off;
axis equal;
axis([-2.1 2.1 -2.1 2.1]);

end

function saveFigure(name)

set(gcf,'PaperPositionMode','auto')
print(name, '-dpng', '-r0');

end

function [c, s1, yy, s2, s1s2] = computeConstantItems(lfsy, eps, angle_y)

% draw C(c, lfs(y))
c = [0 -lfsy]; % center of lfs(y) tangent to y, s1
s1 = [c(1) c(2) + lfsy];

% compute and draw y
syms x y;
% y = circle (s1, eps) intersect circle (c, lfsy)
[x, y] = solve((x-s1(1))^2+(y-s1(2))^2==eps^2, (x-c(1))^2+(y-c(2))^2==lfsy^2);
if x(1)>x(2)
    yy=[double(x(1)) double(y(1))];
else
    yy=[double(x(2)) double(y(2))];
end
alpha = angle_y*atan2(yy(2), yy(1));
yy=[eps*cos(alpha) eps*sin(alpha)];

% compute and draw s2
syms x y;
% s2 = lfsy/2 distance further away on tangent at y (worst case for s1s2)
s2=[yy(1)+eps*cos(2*alpha) yy(2)+eps*sin(2*alpha)];

% draw circle Rs1s2(s1, |s1,s2|)
s1s2=norm(s1 - s2);

end

function [xx, z, startH, endH, posH] = computeVariableItems(lfsy, eps, ratio, angle_x, c, s1, s1s2)

% compute and draw x - on C
syms x y;
% x = circle (s1, ratio*eps) intersect circle (c, 1)
[x, y] = solve((x-s1(1))^2+(y-s1(2))^2==(ratio*eps)^2, (x-c(1))^2+(y-c(2))^2==lfsy^2);

if x(1)<x(2)
    xx=[double(x(1)) double(y(1))];
else
    xx=[double(x(2)) double(y(2))];
end
% x = x-axis rotated ccw by angle with distance ratio from origin
alpha = angle_x*atan2(-xx(2), -xx(1));
xx=[-ratio*eps*cos(alpha) -ratio*eps*sin(alpha)];

% compute H
vecH = [-sin(2*alpha) cos(2*alpha)];
endH = xx + 2*vecH;
startH = xx - 2*vecH;
posH = xx + vecH;

% compute z
syms x y t;
% z = circle (s1, s1s2) intersect line H=c+t*vecH
S = solve(x^2+y^2==s1s2^2, x==xx(1)+t*vecH(1), y==xx(2)+t*vecH(2));
if S.y(1)>S.y(2)
    z = [double(S.x(1)) double(S.y(1))];
else
    z = [double(S.x(2)) double(S.y(2))];
end

end

function drawItems(lfsy, ratio, c, xx, s1, yy, s2, z, s1s2, startH, endH, posH)

% circle Ry
circle(yy(1), yy(2), lfsy, 'black');

% circle Ty
circle(c(1), c(2), lfsy, 'black');

% circle Rx(x, ratio)
circle(xx(1), xx(2), ratio, 'black');

% lfs(y) distance
line([yy(1) yy(1)], [yy(2) yy(2)-lfsy], 'Color', 'black');
text(yy(1), yy(2)-lfsy/2, 'lfs(y)');

% draw edges
line([yy(1) s2(1)], [yy(2) s2(2)], 'Color', 'black');
line([-2 2], [0 0], 'Color', 'black');
line([0 0], [-2 2], 'Color', 'black');
line([xx(1) s1(1)], [xx(2) s1(2)], 'Color', 'black');
plot([startH(1) endH(1)], [startH(2) endH(2)], 'Color', 'red');

% draw points and labels
dot(yy(1), yy(2), 'y', 'black');
dot(s1(1), s1(2), 's1', 'black');
dot(c(1), c(2), 'c', 'black');
dot(s2(1), s2(2), 's2', 'black');
dot(xx(1), xx(2), 'x', 'black');
dot(z(1), z(2), 'z', 'black');

% draw labels
text(posH(1), posH(2), 'H');
text(0, 1.9, 'x');
text(1.9, 0, 'y');
text(s1(1), s1(2)+s1s2, 'Rs1s2');
text(yy(1), yy(2)+1, 'Ry');
text(xx(1), xx(2)+ratio, 'Rx');

end

function drawLfsDisks(lfsy, p, c, ratio, count, str, phase)

max_a = acos(p(2) - c(2));
radius = ratio;

for i=0:count
    if (p(1) > 0)
        a = i/count*(max_a);
    else
        a = -i/count*(max_a);
    end
    syms x y s;
    t = [sin(a) cos(a)];
    t = c + (lfsy + radius)*t;
    if phase == 0
        disk(t(1), t(2), radius, [0.5 0.5 1]);
    else
        circle(t(1), t(2), radius, [0 0 1]);
        text(t(1), t(2) + radius, str, 'Color', 'blue');
        dot(t(1), t(2), '', 'red');
    end
end

end

function drawRtDisks(lfsy, ratio, c, xx, p, str, stepsize, onaxis, phase)

% draw circles Rtx for angles in stepsize steps
max_a = asin(xx(1));

for i=0:stepsize
    a = i/stepsize*max_a;
    syms x y s;
    % point t = circle (p, ratio) intersect line c + s*(0, 1) rotated by angle a
    S = solve(x==s*sin(a),y==-1+s*cos(a),(x-p(1))^2+(y-p(2))^2==ratio^2);
    if double(S.y(1)) > double(S.y(2))
        t=[double(S.x(1)) double(S.y(1))];
    else
        t=[double(S.x(2)) double(S.y(2))];
    end
    if (onaxis == 0) || (t(2)>0) % do not draw circles below x-axis if rolling on it
        if (onaxis == 0)
            radius = norm(t - c) - lfsy; % rolling on lfs(y)
        else
            radius = t(2); % rolling on x-axis (worst case)
        end
        if (phase == 0)
            disktrans(t(1), t(2), radius, [0.5 0.5 1]);
        else
            circle(t(1), t(2), radius, [0 0 1]);
            text(t(1), t(2) + radius, str, 'Color', 'blue');
            dot(t(1), t(2), '', 'red');
        end
    end
end

end

function h = circle(x, y, r, color)

th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit, 'Color', color);

end

function h = disk(x, y, r, color)

th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = fill(xunit, yunit, color);

end

function h = disktrans(x, y, r, color)

th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = fill(xunit, yunit, color, 'FaceAlpha', 0.2);

end

function h = dot(x, y, t, color)

th = 0:pi/50:2*pi;
r=0.02;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = fill(xunit, yunit, color);
text(x + 2*r, y + 2*r, t, 'Color', color);

end

