% read the fixed image
F = readim('fixed.png');
F = im2mat(F,'double');

% read the moved image
G = circshift(F,50,1);
%G = readim('moving.png');
%G = im2mat(G,'double');

% show difference
figCorr=dipshow(G - F);
figCorr.Name='difference';
figCorr.NumberTitle='off';

h_x = 0;
h_y = 0;

[h_x, h_y] = LucasKanade(F, G, h_x, h_y);

F_derivative_x = dx(F);
F_derivative_y = dy(F);

F = F + h_x * F_derivative_x;
F = F + h_y * F_derivative_y;

% show corrected difference
figCorr=dipshow(G - F);
figCorr.Name='corrected difference';
figCorr.NumberTitle='off';


% show corrected F
figCorr=dipshow(F);
figCorr.Name='corrected F';
figCorr.NumberTitle='off';

function [h_x, h_y] = LucasKanade(F, G, h_x, h_y)

    for a = 1:400
        % show fixed image
        %figF=dipshow(F);
        %figF.Name='fixed image';
        %figF.NumberTitle='off';

        % show moving image
        %figG=dipshow(G);
        %figG.Name='moving image';
        %figG.NumberTitle='off';

        % show difference G - F
        %figD=dipshow(G - F);
        %figD.Name='difference image';
        %figD.NumberTitle='off';

        % Squared Difference
        E = sum((G - F).^2);

        % Gradient of F with respect to x
        F_derivative_x = dx(F);
        
        % Gradient of F with respect to y
        F_derivative_y = dy(F);

        % show derivative
        %figFdx=dipshow(F_derivative_x);
        %figFdx.Name='derivative image dx';
        %figFdx.NumberTitle='off';
        
        % show derivative
        %figFdy=dipshow(F_derivative_y);
        %figFdy.Name='derivative image dy';
        %figFdy.NumberTitle='off';
        
        new_hx = sum(transpose(F_derivative_x)* (G - F)) / sum(transpose(F_derivative_x) * F_derivative_x);
        h_x = h_x + new_hx;        

        new_hy = sum(transpose(F_derivative_y) * (G - F)) / sum(transpose(F_derivative_y) * F_derivative_y);
        h_y = h_y + new_hy;

        F = F + h_x * F_derivative_x;
        F = F + h_y * F_derivative_y;

        % show corrected F
        %figCorr=dipshow(G - F);
        %figCorr.Name='corrected difference';
        %figCorr.NumberTitle='off';

        %close all
    end
    
    
end