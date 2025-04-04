Image = imread ('blurGirl.png');
OriginalImage = double (Image)/255;
g = B(OriginalImage)+0.01*randn(size(OriginalImage));
u1 = DeBlur_TV2(g);
subplot(1,2,1), imshow(g)
title('Blured Image');
subplot(1,2,2), imshow(u1)
title('DeBlured Using TV2 Image');

%De-Blur Using TV
function [u]=DeBlur_TV2(g)
    Sigma = 5/37;
    Tau = 5/37;
    lambda = 0.002;
    [m, n] = size(g);
    u = rand(m,n);
    uBar = rand(m,n);
    y1 = rand(m,n);
    y2 = rand(m,n,2);
    Bu = B(u);
    Gradu = grad(u);
    e = 0.01;
    PD = CondTV2(Bu,g,Gradu,lambda,y1);
    while abs(PD)>=e
        disp(abs(PD) - e);
        uHat = u;
        y1 = (y1 + Sigma * (B(uBar) - g)) / (1 + Sigma);
        y2 = POmega1(y2 + Sigma * grad(uBar), lambda);
        u = uHat - Tau * (BStar(y1) - div(y2));
        uBar = 2*u - uHat;
        Bu = B(u);
        Gradu = grad(u);
        PD = CondTV2(Bu,g,Gradu,lambda,y1);
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
function [x] = grad(A)
    [m,n] = size(A);
    x = zeros(m,n,2);
    x( : , : , 1 ) = DxPlus(A);
    x( : , : , 2 ) = DyPlus(A);
end

%Divergence Function
function [x] = div(A)
    x = DxMinus(A( : , : , 1 )) + DyMinus(A( : , : , 2 ));
end

%ProxFstar Function
function [x] = POmega1(P,Sigma)
    [m,n,~] = size(P);
    Q = max( Sigma*ones(m,n) , sqrt(P( : , : , 1 ).^2 + P( : , : ,2).^2));
    x = Sigma*P./Q;
end

%Matrix Norm-2
function [x] = M2(P)
    x = sqrt(sum(sum((P.^2))));
end

%Norm 1,2 Computation Function
function [x] = M12(P)
	x = sum(sum(sqrt((P( : , : , 1 ).^2) + (P( : , : ,2).^2))));
end

%Convolution
function [x] = B(u)
    % Gaussian Blur Kernel:
    ker=imread('motion_kernel_3.tif');
    ker=double(ker);
    x = conv2(u,ker,'same');
end

%ConvolutionStar
function [x] = BStar(v)
    % Gaussian Blur Kernel:
    ker=imread('motion_kernel_3.tif');
    ker=double(ker);
    x = conv2(v,ker(end:-1:1,end:-1:1),'same');
end

% %Condition Check Function
function [x] = CondTV2(Bu,g,Gradu,lambda,y1)
    x = (0.5)*((M2(Bu-g))^2) + lambda*M12(Gradu) + sum(sum(y1.*g)) + (0.5)*(M2(y1)^2);
end