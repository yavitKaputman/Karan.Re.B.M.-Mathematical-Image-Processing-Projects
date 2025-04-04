Image = imread ('noisy01.jpg');
g1 = double (Image)/255;
Image = rgb2gray(imread ('noisy02.jpg'));
g2 = double (Image)/255;
Image = imread ('noisy03.png');
g3 = double (Image)/255;
u1 = TGV(g1);
u2 = TGV(g2);
u3 = TGV(g3);
subplot(3,2,1), imshow(g1)
title('Input Image');
subplot(3,2,2), 
imshow(u1)
title('Output Image');
subplot(3,2,3), imshow(g2)
title('Input Image');
subplot(3,2,4), imshow(u2)
title('Output Image');
subplot(3,2,5), imshow(g3)
title('Input Image');
subplot(3,2,6), imshow(u3)
title('Output Image');

%Den-TGV Function
function [u]=TGV(g)
    Sigma = 5/37;
    Tau = 5/37;
    alpha0 = 0.1;
    alpha1 = 0.1;
    [m, n] = size(g);
    u = zeros(m,n);
    w = zeros(m,n,2);
    uBar = zeros(m,n);
    wBar = zeros(m,n,2);
    p = zeros(m,n,2);
    q = zeros(m,n,3);
    Gradu = Grad(u);
    Divp = Div(p);
    Ew = E(w);
    i = 1;
    Epsilon = 0.2;
    Cond = Condition(u,w,g,Ew,Gradu,Divp,alpha0,alpha1);
    while Cond >= Epsilon
        disp(Cond - Epsilon);
        uHat = u;
        wHat = w;
        p = POmega1(p+Sigma*(Grad(uBar)-wBar), alpha1);
        q = POmega2(q+Sigma*(E(wBar)), alpha0);
        u = ProxG1(uHat,g,Divp,Tau);
        w = ProxG2(wHat,p,divPrime(q),Tau);
        uBar = 2*u - uHat;
        wBar = 2*w - wHat;
        Gradu = Grad(u);
        Divp = Div(p);
        Ew = E(w);
        if i == 2000
            break
        else
            i = i + 1;
        end
        Cond = Condition(u,w,g,Ew,Gradu,Divp,alpha0,alpha1);
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
function [x] = Grad(A)
    [m,n] = size(A);
    x = zeros(m,n,2);
    x( : , : , 1 ) = DxPlus(A);
    x( : , : , 2 ) = DyPlus(A);
end

%Divergence Function
function x = Div(A)
    x = DxMinus(A( : , : , 1 )) + DyMinus(A( : , : , 2 ));
end

%Epsilon Function
function [v] = E(w)
[m,n,~]=size(w);
v = zeros(m,n,3);
v(:,:,1) = DxPlus(w(:,:,1));
v(:,:,2) = (0.5)*(DyPlus(w(:,:,1)) + DxPlus(w(:,:,2)));
v(:,:,3) = DyPlus(w(:,:,2));
end

%divPrime Function
function [v] = divPrime(y)
[m,n]=size(y(:,:,1));
v = zeros(m,n,2);
v(:,:,1) = DxMinus(y(:,:,1))+DyMinus(y(:,:,2));
v(:,:,2) = DxMinus(y(:,:,2))+DyMinus(y(:,:,3));
end

%POmega Functions
function [x] = POmega1(p,Sigma)
    [m,n,~] = size(p);
    Q = max( Sigma*ones(m,n) , sqrt((p( : , : , 1 )).^2 + (p( : , : ,2)).^2));
    x = Sigma*p./Q;
end
function [x] = POmega2(q,Sigma)
    [m,n,~] = size(q);
    Q = max( Sigma*ones(m,n) , sqrt((q( : , : , 1 )).^2 + 2*(q( : , : ,2)).^2 + (q( : , : ,3)).^2));
    x = Sigma*q./Q;
end

%ProxG
function [u] = ProxG1(uHat,g,Divp,Tau)
u = (uHat + Tau*(Divp+g))/(1+Tau);
end
function [u] = ProxG2(wHat,p,divPrimeq,Tau)
u = (wHat + Tau*(divPrimeq+p));
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
function [x] = Condition(u,w,g,Ew,Gradu,Divp,alpha0,alpha1)
x = (0.5)*((M2(u-g))^2) + alpha1*TV(Gradu - w) + alpha0*TV(Ew) + sum(sum(Divp.*g)) + (0.5)*M2(Divp)^2;
end