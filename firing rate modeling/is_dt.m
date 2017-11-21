function [ iss, ds ] = is_dt( is_0, w, u, tau_s )
%IS_DT I_s(t) calculation w/ dI_s/dt
%   is_0: first I_s
%   w: weights, N_u X # postsynaptic neurons
%   u: presynaptic firing rates wrt time
%   tau_s: time constant
%   ds: deltas
%   iss: I_s's

ds = zeros(size(u));
iss = zeros(size(u));
is = is_0;
for i = 1 : size(u,2)
    ds(i) = dis(is, w, u(i), tau_s);
    is = is + ds(i);
    iss(i) = is;
end

end