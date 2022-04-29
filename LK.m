img = zeros(256,256);
img(25:50, 25:50) = 255;
in1 = dipshow(img);
in1.Name='orig';

img2 = zeros(256,256);
img2(30:55, 30:55) = 255;
in2 = dipshow(img2);
in2.Name='shift';

h = [0; 0];

h = LucasKanade(img, img2, h, img2)

regist = imtranslate(img2,(h)');
in3 = dipshow(regist);
in3.Name='regist';

function h = LucasKanade(img, img2, h, extra)
    
    for iter = 1:100
        
        imgDx = dx(img2, 1);
        imgDy = dy(img2, 1);
        
        J_prime_xx = imgDx.*imgDx;
        J_prime_xy = imgDx.*imgDy; 
        J_prime_yy = imgDy.*imgDy;

        H = [sum(J_prime_xx) sum(J_prime_xy); sum(J_prime_xy) sum(J_prime_yy)];
        
        img2 = imtranslate(img2,h');
        diff = img - img2;
                
        if abs(sum(diff)) < 0.5
            break
        end
        num = [-sum(imgDx * diff) -sum(imgDy * diff)];

        deltaH = num / H;

        h = h + transpose(deltaH);
    end
    disp(iter)
end