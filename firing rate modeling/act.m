function y = act( is, delta )
%F F(I_s) activation function
%   is: I_s
%   delta: threshold

if nargin < 2
    delta = .5;
end

y = halfwave_rect(is - delta);
% y = sigmoid(is);

end