function imgArray = xydata2arr(X, Y, data)
% XYDATA2ARR   Convert an XY labeled dataset into an X by Y array
% X - integer col vector of all X indices  (X1; X2; … Xi; … Xn)
% Y - integer col vector of all Y indices  (Y1; Y2; … Yi; … Yn)
% data - col vector of XY data  (f(X1,Y1); f(X2,Y2); … f(Xi,Yi); … f(Xn, Yn))

imgArray=accumarray([Y X],data);
