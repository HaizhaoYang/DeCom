function y = gen_shape(x,num)
% Generate wave shapes

x = x - floor(x);
y = zeros(size(x));
if num == 1
    loc = find(x<0.5);
    y(loc) = sqrt(0.25^2-(x(loc)-0.25).^2);
    loc = find(x>=0.5);
    y(loc) = -sqrt(0.25^2-(x(loc)-0.75).^2);
    y = y - mean(y);
    y = y/sqrt(sum(abs(y).^2)/length(y));
else if num == 2
    loc = find(x<=0.25);
    y(loc) = x(loc);
    loc = find(x>0.25 & x<=0.75);
    y(loc) = 0.5-x(loc);
    loc = find(x>0.75);
    y(loc) = x(loc)-1;
    y = y - mean(y);
    y = y/sqrt(sum(abs(y).^2)/length(y));
    else if num ==3
        loc = find(x<=0.5);
        y(loc) = 1;
        loc = find(x>0.5);
        y(loc) = -1;
        y = y - mean(y);
        y = y/sqrt(sum(abs(y).^2)/length(y));
        else if num == 4
            loc = find(x<=0.25);
            y(loc) = 1;
            loc = find(x>0.25 & x<=0.5);
            y(loc) = 3-8*x(loc);
            loc = find(x>0.5 & x<= 0.75);
            y(loc) = -1;
            loc = find(x>0.75);
            y(loc) = 8*x(loc)-7;
            y = y - mean(y);
            y = y/sqrt(sum(abs(y).^2)/length(y));
            else if num ==5
                loc = find(x<=0.45);
                y(loc) = x(loc);
                loc = find(x>0.45 & x<=0.55);
                y(loc) = 4.5-9*x(loc);
                loc = find(x>0.55);
                y(loc) = x(loc)-1;
                y = y - mean(y);
                y = y/sqrt(sum(abs(y).^2)/length(y));
                else if num == 6
                    x = 2*pi*x;
                    y = (cos(0.8*cos(x))-tan(x).*sin(0.8*cos(x))).*cos(x);
                    y = y - mean(y);
                    y = y/sqrt(sum(abs(y).^2)/length(y));
                    else
                        loc = find(x<=0.3);
                        y(loc) = 0.5;
                        loc = find(x>0.3 &x<=0.6);
                        y(loc) = -0.2;
                        loc = find(x>0.6);
                        y(loc) = 0.2;
                        y = y - mean(y);
                        y = y/sqrt(sum(abs(y).^2)/length(y));
                    end
                end
            end
        end
    end
end