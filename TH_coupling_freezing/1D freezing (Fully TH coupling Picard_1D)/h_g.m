function [ h_g ] = h_g(H, L)
% Calculate the gradient of water head
h_g = (H(1) - H(2))/L ;

end

