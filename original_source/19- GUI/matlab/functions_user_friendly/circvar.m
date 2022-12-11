function[cv]=circvar(orients, resps)
%orients = orientations in radians
% resps = resps (in order of orients) 
%orients=orients*(2*pi/180);
exponent=exp(1i*2*orients); % create normal vectors with angle orients
num=resps*exponent'; % weight the vectors by response size
denom=sum(resps); %sum the total responses
R=num./denom; %take ratio
cv=1-abs(R); % one minus because literature is always this direction ( so 1 = more variable and 0 = more selective)