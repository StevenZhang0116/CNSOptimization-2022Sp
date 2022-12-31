clear all;close all;clc
format long;
%% (5.114) for data set 1
one = ones(10,1);
load('hw4data1.mat','W');
cvx_begin sdp
    variable v1(10);
    dual variable x;
    maximize (-(one')*v1);
    x: W+diag(v1)>=0;
cvx_end
%% (5.114) for data set 2
one = ones(50,1);
load('hw4data2.mat','W');
cvx_begin sdp
    variable v2(50);
    dual variable lambda;
    maximize (-(one')*v2);
    lambda: W+diag(v2)>=0;
cvx_end
%% (5.115) for data set 1
one = ones(10,1);
load('hw4data1.mat','W');
cvx_begin sdp
    variable X1(10,10)
    dual variables y1 z1;
    minimize trace(W*X1);
    subject to 
               y1: X1>=0;
               z1: diag(X1) == one;
cvx_end
rankofprimal1 = compute_rank(X1);
rankofdual1 = compute_rank(y1);

%% (5.115) for data set 2
one = ones(50,1);
load('hw4data2.mat','W');
cvx_begin sdp
    variable X2(50,50)
    dual variables y2 z2;
    minimize trace(W*X2);
    subject to 
               y2: X2>=0;
               z2: diag(X2) == one;
cvx_end
rankofprimal2 = compute_rank(X2);
rankofdual2 = compute_rank(y2);

% display(rankofprimal1)
% display(rankofdual1)
% display(rankofprimal2)
% display(rankofdual2)

%% (c) carrying out the algorithm on p.1120

%% data set 1
load('hw4data1.mat','W');
cvx_begin sdp
    variable X1(10,10)
    dual variables y1 z1;
    minimize trace(W*X1);
    subject to 
               y1: X1>=0;
               z1: diag(X1) == ones(10,1);
cvx_end
V1 =chol(X1);
r = ones(10,1)*(1/(10^(1/2)));
% partition
p1 = zeros(10,1);
for k = 1:10
    if(r'*V1(:,k)>=0)
        p1(k) = 1;
    else
        p1(k) = -1;
    end
end
% sum
S1 = 0;
for j = 2:10
    for i=1:j-1
        S1 = S1 + 1/2*(W(i,j)-W(i,j)*(p1(i)*p1(j)));
    end
end
display(S1);

%% data set 2
load('hw4data2.mat','W');
cvx_begin sdp
    variable X2(50,50)
    dual variables y2 z2;
    minimize trace(W*X2);
    subject to 
               y2: X2>=0;
               z2: diag(X2) == ones(50,1);
cvx_end
V2 = chol(X2);
% partition
r = ones(50,1)*(1/(50^(1/2)));
p2 = zeros(50,1);
for k = 1:50
    if(r'*V2(:,k)>=0)
        p2(k) = 1;
    else
        p2(k) = -1;
    end
end

% sum
S2 = 0;
for j = 2:50
    for i=1:j-1
        S2 = S2 + 1/2*(W(i,j)-W(i,j)*(p2(i)*p2(j)));
    end
end
display(S2);

%% brute force method
load('hw4data1.mat','W');
OptimalforSDP = realmax;
MaxCut = realmin;
x1 = ones(10,1);
for i1= 1:2
    x1(1) = (-1)^i1;
for i2 = 1:2
    x1(2) = (-1)^i2;
for i3 = 1:2
    x1(3) = (-1)^i3;
for i4 = 1:2
    x1(4) = (-1)^i4;
for i5 = 1:2
    x1(5) = (-1)^i5;
for i6 = 1:2
    x1(6) = (-1)^i6;
for i7 = 1:2
    x1(7) = (-1)^i7;  
for i8 = 1:2
    x1(8) = (-1)^i8;    
for i9 = 1:2
    x1(9) = (-1)^i9;
for i10 = 1:2
    x1(10) = (-1)^i10;
    if(x1'*W*x1)<OptimalforSDP
        OptimalforSDP = x1'*W*x1;
    end
sum = 0;
for j = 2:10
    for i=1:j-1
        sum = sum + (1/2)*(W(i,j)*(1 - x1(i)*x1(j)));
    end
end
  if sum >= MaxCut
      MaxCut = sum;
  end
end
end
end
end
end
end
end
end
end
end

display(MaxCut);
display(OptimalforSDP);

%% compute rank of a matrix
function r = compute_rank(A)
    e = eig(A);
    count = 0;
    for k = 1:length(e)
        if e(k) < 1e-5
            count = count + 1;
        end   
    end
r = length(A) - count;
end


    
    
