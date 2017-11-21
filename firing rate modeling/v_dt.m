function [ vs, ds ] = v_dt( v_0, is_t, tau_r )
%V_DT v calculation w/ dv/dt (postsyn firing rate)
%   v_0: first v
%   is_t: I_s wrt time
%   tau_r: time constant

vs = zeros(size(is_t));
ds = zeros(size(is_t));
v = v_0;
for i = 1 : size(is_t,2)
    ds(i) = dv(v, is_t(i), tau_r);
    v = v + ds(i);
    vs(i) = v;
end

end