function b = skew(a)
%Skew calculates the 3x3 skew symetric matrix of a 3x1 vector

if size(a) ~= [3 1]
    error('invalid vector dimesions')
else
    b = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end 

end