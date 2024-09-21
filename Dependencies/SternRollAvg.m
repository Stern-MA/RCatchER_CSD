%this function take a rolling average of an input vector (vec1)
%it does this by taking the average of m elements before and after
%a given element, thus smoothing out the input vector
%Note: the bookends of the vector (m number of elements) are constrained to
%what can be averaged only (i.e. element 1 will only be averaged with m
%elements following it as there are none preceding it, and the last element
%will only be averaged with m elements preceding it as there are none
%following it
function [RA]= SternRollAvg(vec1,m)
if size(vec1,1)==1
    tempVec1=[nan(1,m),vec1,nan(1,m)];
elseif size(vec1,2)==1
    tempVec1=[nan(1,m),vec1',nan(1,m)];
else
    disp('Error, please enter vector')
end
RA=nan(length(vec1),1);
for ii = (m+1):(m+length(vec1))
    RA(ii-m)=nanmean(tempVec1((ii-m):(ii+m)));
end
clear tempVec1
clear ii