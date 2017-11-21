% constants / initial values
is_0 = 0;
v_0 = 0;
w = 5;
tau_s = 2;
tau_r = 30;
% inputs
u = [zeros(1,10) ones(1,10) zeros(1,1)...
    ones(1,20) zeros(1,10) ones(1,4) zeros(1,30)];

% simulate I_s
iss = is_dt(is_0, w, u, tau_s);
subplot(1,2,1)
plot(iss, 'LineWidth', 2)
hold on
ax = gca;
plot(find(u>0), ax.YLim(2) * ones(size(find(u>0))),...
    '.')
hold off; ph.prefs

% simulate v
vs = v_dt(v_0, iss, tau_r);
subplot(1,2,2)
plot(vs, 'LineWidth', 2)
hold on
ax = gca;
plot(find(u>0), ax.YLim(2) * ones(size(find(u>0))),...
    '.')
hold off; ph.prefs
