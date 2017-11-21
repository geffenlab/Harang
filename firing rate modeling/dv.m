function d = dv( v, is, tau )
%DV delta v (firing rate)

d = (-1*v + act(is)) / tau;

end