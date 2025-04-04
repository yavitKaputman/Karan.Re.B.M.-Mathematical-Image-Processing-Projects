Image = rgb2gray(imread('TestSubject.jpg'));
g = double (Image)/255;
[m , n] = size(g);
if rem(m,2) ~= 0
    g = [g;g(end,:)];
end
if rem(n,2) ~= 0
    g = [g g(:,end)];
end
u = DS(g);
v = DSstar(u);
subplot(1,3,1), imshow(g)
title('Input Image');
subplot(1,3,2), imshow(u)
title('Ds Image');
subplot(1,3,3), imshow(v)
title('DSstar Image');
u = rand(2000,2000);
v = rand(1000,1000);
e = sum(sum(DS(u).*v)) - sum(sum(u.*DSstar(v)));
disp("Error Is : " + e);

function u = DS(g)
u = 1/4*(g(1:2:end,1:2:end) + g(2:2:end,1:2:end) + g(1:2:end,2:2:end) + g(2:2:end,2:2:end));
end

function g = DSstar(u)
[m,n] = size(u);
g = zeros(2*m,2*n);
g(1:2:end,1:2:end) = 1/4*u;
g(2:2:end,1:2:end) = 1/4*u;
g(1:2:end,2:2:end) = 1/4*u;
g(2:2:end,2:2:end) = 1/4*u;
end