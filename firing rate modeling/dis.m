function d = dis( is, w, u, tau )
%DIS delta I_s

d = (-1*is + w' * u) / tau;

end