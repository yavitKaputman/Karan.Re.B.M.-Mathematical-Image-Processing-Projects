Image = imread ('noisy01.jpg');
g1 = double (Image)/255;
Image = rgb2gray(imread ('noisy02.jpg'));
g2 = double (Image)/255;
Image = imread ('noisy03.png');
g3 = double (Image)/255;
u1 = ROF(g1);
u2 = ROF(g2);
u3 = ROF(g3);
subplot(3,2,1), imshow(g1)
title('Input Image');
subplot(3,2,2), imshow(u1)
title('Output Image');
subplot(3,2,3), imshow(g2)
title('Input Image');
subplot(3,2,4), imshow(u2)
title('Output Image');
subplot(3,2,5), imshow(g3)
title('Input Image');
subplot(3,2,6), imshow(u3)
title('Output Image');

%ROF Function
function [u]=ROF(g)
    Sigma = 5/37;
    Tau = 5/37;
    lambda = 0.1;
    [m, n] = size(g);
    u = rand(m,n);
    uBar = rand(m,n);
    p = rand(m,n,2);
    Gradu = Gradient(u);
    Divp = Divergence(p);
    i = 1;
    Epsilon = 0.1;
    Cond = Condition(u,g,Divp,Gradu,lambda);
    while Cond >= Epsilon
        disp(Cond - Epsilon);
        uHat = u;
        p = ProxFstar(p+Sigma*Gradient(uBar),Sigma);
        u = ProxG(uHat,Divp,g,Tau);
        uBar = 2*u - uHat;
        Gradu = Gradient(u);
        Divp = Divergence(p);
        if i == 2000
            break
        else
            i = i + 1;
        end
        Cond = Condition(u,g,Divp,Gradu,lambda);
    end
end

% DxPlus Function
function [Dx] = DxPlus(D)
    [m,n] = size(D);
    Dx = zeros(m,n);
    Dx( :  , 1 : n - 1 ) = D( : , 2 : n ) - D( : , 1 : n-1 );
end

% DyPlus Function
function [Dy] = DyPlus(D)
    Dy = DxPlus(D')';
end

% DxMinus Function
function [Dx] = DxMinus(D)
    [m,n] = size(D);
    Dx = zeros(m,n);
    Dx( : , 1 ) = D( : , 1);
    Dx( : , 2 : n - 1 ) = D( : , 2 : n - 1 ) - D( : , 1 : n - 2 );
    Dx( : , n ) = -D( : , n - 1 );
end

% DyMinus Function
function [Dy] = DyMinus(D)
    Dy = DxMinus(D')';
end

% Gradient Function
function GradA = Gradient(A)
    [m,n] = size(A);
    GradA = zeros(m,n,2);
    GradA( : , : , 1 ) = DxPlus(A);
    GradA( : , : , 2 ) = DyPlus(A);
end

%Divergence Function
function DivA = Divergence(A)
    DivA = DxMinus(A( : , : , 1 )) + DyMinus(A( : , : , 2 ));
end

%ProxFstar Function
function [POmega] = ProxFstar(P,Sigma)
    [m,n,~] = size(P);
    Q = max( Sigma*ones(m,n) , sqrt((P( : , : , 1 )).^2 + (P( : , : ,2)).^2));
    POmega = Sigma*P./Q;
end

%ProxG
function [u] = ProxG(uHat,Divp,g,Tau)
u = (uHat + Tau*(Divp+g))/(1+Tau);
end

%Matrix Norm-2
function [x] = M2(P)
    x = sqrt(sum(sum((P.^2))));
end

%TV2 Computation Function
function [x] = TV(P)
    x = sum(sum(sqrt((P( : , : , 1 ).^2) + (P( : , : ,2).^2))));
end

%Condition Check Function
function [x] = Condition(u,g,Divp,Gradu,lambda)
x = (0.5)*((M2(u-g))^2) + lambda*TV(Gradu) + sum(sum(Divp.*g)) + (0.5)*M2(Divp)^2;
end