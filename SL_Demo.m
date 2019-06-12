function [] = SL_Demo(x,y,q,sl)
    
    %constants
    phi = 0.2;
    h = 20;
    u_w = 1;
    u_o = 10;
    k = 150;
    r_w = 0.25;

    %for conversion
    fac = 1/(2*pi*h*phi);
    
    %invert sign of injectgion (for convinience)
    q = -1*q;
    
    %define dt and number of SL's
    dt = 1; %2/log2(sl);
    %q = q/sl;
    
    %velocity to flux convertion factor
    v2q = 2*pi*h*phi*r_w/sl;

    %location of streamlines around the well
    pp = [-pi:2*pi/sl:pi];
    xx = cos(pp);
    yy = sin(pp);

    %Trace SL's from injection well
    i = 1;
    v_hist = 0;
    conv_sl = 0;
    err_SL_bundle = zeros(size(x,2)-1,1);
    for k = 1:sl
        xp = xx(k);
        yp = yy(k);
        pt = 1;
        SL(i,pt,1) = xp;
        SL(i,pt,2) = yp;
        reached_prod_well = 0;
        clear v_hist;
        while(reached_prod_well == 0 && pt < 20000*log2(sl))
            %calc vx and vy
            vx = 0;
            vy = 0;
            for j = 1:size(x,2)
                denom = ((xp - x(j))^2 + (yp - y(j))^2);
                vx = vx + q(j)*(xp - x(j))/denom;
                vy = vy + q(j)*(yp - y(j))/denom;
            end
            v_hist(pt) = sqrt(vx^2+vy^2);
            vx = fac*vx;
            vy = fac*vy;
            
            xp = xp + vx*dt;
            yp = yp + vy*dt;
            
            pt = pt + 1;
            SL(i,pt,1) = xp;
            SL(i,pt,2) = yp;
            reached_prod_well = isCloseToProd(xp, yp, x, y);
        end
        pt
        if(reached_prod_well > 1)
                i
                conv_sl = conv_sl + 1
                err_SL(conv_sl) = v2q*(v_hist(end)-v_hist(1));
                err_SL_bundle(reached_prod_well-1) = err_SL_bundle(reached_prod_well-1) + err_SL(conv_sl);
                disp('prod well reached');
        else
                i
                disp('prod well Not reached');
        end

        i = i+1;
    end
    
    %Trace SL's from production wells
    q = -1*q;
    for k = 2:size(x,2)
        conv_sl = 0;
        for i = 1:sl
            xp = xx(i)+x(k);
            yp = yy(i)+y(k);
            pt = 1;
            SL(i+(k-1)*sl,pt,1) = xp;
            SL(i+(k-1)*sl,pt,2) = yp;
            reached_inj_well = 0;
            clear v_hist;
            while(reached_inj_well == 0 && pt < 20000*log2(sl))
                %calc vx and vy
                vx = 0;
                vy = 0;
                for j = 1:size(x,2)
                    denom = ((xp - x(j))^2 + (yp - y(j))^2);
                    vx = vx + q(j)*(xp - x(j))/denom;
                    vy = vy + q(j)*(yp - y(j))/denom;
                end
                v_hist(pt) = sqrt(vx^2+vy^2);
                
                vx = fac*vx;
                vy = fac*vy;
                
                xp = xp + vx*dt;
                yp = yp + vy*dt;
                
                pt = pt + 1;
                SL(i+(k-1)*sl,pt,1) = xp;
                SL(i+(k-1)*sl,pt,2) = yp;
                reached_inj_well = isClose(xp, yp, x, y);
            end
            pt
            if(reached_inj_well == 1)
                i+(k-1)*sl
                conv_sl = conv_sl + 1
                err_SL(conv_sl) = v2q*(v_hist(end)-v_hist(1));
                err_SL_bundle(k-1) = err_SL_bundle(k-1) - err_SL(conv_sl);
                disp('inj well reached');
            else
                i+(k-1)*sl
                disp('inj well Not reached');
            end
        end
    end
    
    %plot SL's from injector
    for i = 1:sl%
        plot(SL(i,2:end-1,1),SL(i,2:end-1,2), '.','MarkerSize',1,'MarkerEdgeColor','b')
        hold on
    end
    
    for i = sl+1:size(x,2)*sl
        plot(SL(i,2:end-1,1),SL(i,2:end-1,2), '.','MarkerSize',1,'MarkerEdgeColor','r')
        hold on
    end
    
    %plot SL's at 5 different times...
    timespan = 1500;
    colors = ['b';'r';'c';'g';'k'];
    figure;
    
    %plot streamlines..
    for j = 1:5
        for i = 1:sl
            front(j,i,1)= SL(i,timespan*j,1);
            front(j,i,2)= SL(i,timespan*j,2);
            hold on
        end
        front(j,i+1,1)= SL(1,timespan*j,1);
        front(j,i+1,2)= SL(1,timespan*j,2);
        plot(front(j,:,1), front(j,:,2), [colors(j),'-'], 'LineWidth',3);
        hold on
    end
    
    legend('1500','3000','4500','6000','7500');
    
    plot(x(2:end),y(2:end),'x','MarkerSize',8,'MarkerEdgeColor','r');
    hold on
    plot(x(1),y(1),'o','MarkerSize',8);
    hold on
    %plot streamlines..
    for j = 1:5
        for i = 1:sl
            plot(SL(i,timespan*(j-1)+1:timespan*j,1),SL(i,timespan*(j-1)+1:timespan*j,2), '.','MarkerSize',1,'MarkerEdgeColor',colors(j));
            front(j,i,1)= SL(i,timespan*j,1);
            front(j,i,2)= SL(i,timespan*j,2);
            hold on
        end
    end
    
    %Error Norms
    figure;
    plot(err_SL,'o');
    figure;
    plot(err_SL_bundle,'o');
    norm(err_SL,inf)
    norm(err_SL_bundle,inf)
    
end

function [close] = isCloseToProd(xp, yp, x, y)
    tol = 0.1;
    close = 0;
    for i = 2:size(x,2)
         if(sqrt((xp-x(i))^2+(yp-y(i))^2) < tol)
             close = i;
         end
    end
end

function [close] = isClose(xp, yp, x, y)
    tol = 0.1;
    close = 0;
    i = 1;
    if(sqrt((xp-x(i))^2+(yp-y(i))^2) < tol)
        close = 1;
    end
end
