function y = halfwave_rect( x )
%HALFWAVE_RECT "F()" half-wave rectification

y = x;
y(y<0) = 0;

end