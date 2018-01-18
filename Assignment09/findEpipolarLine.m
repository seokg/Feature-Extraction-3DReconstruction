function [ slope,b ] = findEpipolarLine( epipoleline )
% Finde the slope and b of the epipolar line

slope = zeros(1,length(epipoleline));
b = zeros(1,length(epipoleline));
for i = 1:length(epipoleline)
slope(1,i) = - epipoleline(1,i) / epipoleline(2,i);
b(1,i) = - epipoleline(3,i) / epipoleline(2,i);
end

end

