function [ hs ] = relief( X,Y,t )
global test
global length h0

if test == 0
    r1=sqrt((X-.5*length).^2+(Y-.5*length).^2);
    hs=.8*h0.*exp(-(r1.^2./(.5*length)));
elseif test == 1
    r1=sqrt((X-.5*length).^2+(Y-.5*length).^2);
    hs=.8*h0.*exp(-(r1.^2./(.5*length)));
else
    zeros(size(X));
end

end

