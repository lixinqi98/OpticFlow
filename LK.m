% img = zeros(256,256);
% img(25:50, 25:50) = 255;
% in1 = dipshow(img);
% in1.Name='orig';
% 
% img2 = zeros(256,256);
% img2(26:51, 30:55) = 255;
% in2 = dipshow(img2);
% in2.Name='shift';

F = readim('fixed.png');
F = im2mat(F,'double');
F = F(:,:,1);

% read the moved image
G = imtranslate(F, [0.5, 0]);

% show the difference
figure; imshowpair(F,G,'ColorChannels','red-cyan');
title('Color composite (Original = red, Shifted = cyan)');
%%
h = [0; 0];

[h, record_h] = LucasKanade(G, F, h);

G_ = imtranslate(F, h');
figure; imshowpair(G,G_,'ColorChannels','red-cyan');
title('Color composite (Original = red, Registered = cyan)');

disp(h')
figure, plot(record_h(:,1)); title('changes in x')
figure, plot(record_h(:,2)); title('changes in y')

% regist = imtranslate(img2,(h)');
% in3 = dipshow(regist);
% in3.Name='regist';

function [h, record] = LucasKanade(img, img2, h)
%     err_0 = immse(img,img2);
    for iter = 1:500
%         TODO: whether sigma affects?
%         err = immse(img,img2);
        
        img2 = imtranslate(img2,h'); 
        
        imgDx = dx(img2, 1);
        imgDy = dy(img2, 1);
        
        J_prime_xx = imgDx.*imgDx;
        J_prime_xy = imgDx.*imgDy; 
        J_prime_yy = imgDy.*imgDy;

        H = [sum(J_prime_xx) sum(J_prime_xy); sum(J_prime_xy) sum(J_prime_yy)];
        
        
        diff = img - img2;     
        
        num = [-sum(imgDx * diff) -sum(imgDy * diff)];
        
        deltaH = num / H;
        
        if isnan(deltaH)
            disp("The H is not inversible")
            break
        end
        
%         if sum(deltaH) < 1e-8
%             disp(deltaH)
%             disp("Convergent")
%             break
%         end
        
        h = h + transpose(deltaH);
        record(iter,:) = h';
    end
    disp(iter)
end