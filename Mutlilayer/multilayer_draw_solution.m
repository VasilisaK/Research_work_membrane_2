function multilayer_draw_solution

%A0 = -0.395617582102597;
%A1 = 0.275304624145246;
%R0 = -1.0;

%A0 = -1.988829934363611e-07;
%A1 = -1.400114748948522e-07;
%A2 = -4.249116508581264e-07;
%R0 = 1.000001;

A0 = -0.081392384473940;
A1 = -0.240888873266162;
A2 = -0.459153316022211;
R0 = 1.77452;

fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Multilayer\res.txt','w');
fprintf(fileID,'\n%f\t\t%f\t\t%f\t',A0,A1,R0);
fclose(fileID);

% solution plot
%draw_solution_field(A0,A1,R0);
draw_solution_field_2(A0,A1,A2,R0);

end

% conformal mapping, new polar coordinates
function [rho, Teta] = calc_cyl_2d(r,Fi_rad)

r0 = 3.95; % domain under cargo radius
Rm = 50; % membrane radius
a = 10; % distance between membrane center and cargo center

Qplus = Rm^2 + a^2 - r0^2 + sqrt( (Rm^2 + a^2 - r0^2)^2 - 4*a^2*Rm^2 );
Qminus = Rm^2 + a^2 - r0^2 - sqrt( (Rm^2 + a^2 - r0^2)^2 - 4*a^2*Rm^2 );
k = Qplus/(2*a*Rm);

rho = zeros(size(r));
Teta = zeros(size(Fi_rad));

for i = 1:size(rho,1)
    for j = 1:size(Teta,2)
        rho_up = 4*a^2*r(i,j)*r(i,j) + 4*a*Qminus*r(i,j)*cos(Fi_rad(i,j)) + Qminus^2;
        rho_down = 4*a^2*r(i,j)*r(i,j) + 4*a*Qplus*r(i,j)*cos(Fi_rad(i,j)) + Qplus^2;
        rho(i,j) = k*sqrt( rho_up/rho_down );
        Teta_up = r(i,j)*sqrt( ( Rm^2+a^2-r0^2 )^2 - 4*Rm^2*a^2 )*sin(Fi_rad(i,j));
        Teta_down = a*(Rm^2+r(i,j)*r(i,j)) + r(i,j)*(Rm^2+a^2-r0^2)*cos(Fi_rad(i,j));
        Teta(i,j) = atan2( Teta_up,Teta_down );
    end
end


end

function draw_solution_field(A0,A1,R0)

R = 1.0;
Nx = 51;
Ny = 51;
xn = linspace(-50,50,Nx);
yn_t = linspace(-50,50,Ny);
yn = yn_t';
coords2d_r = zeros(size(yn,1),size(xn,2));
coords2d_fi = zeros(size(yn,1),size(xn,2));


for i = 1:size(yn,1)
    for j = 1:size(xn,2)
        coords2d_fi(i,j) = atan2(yn(i,1),xn(1,j));
    end
end
 
for i = 1:size(yn,1)
    for j = 1:size(xn,2)
        coords2d_r(i,j) = sqrt(xn(1,j)^2 + yn(i,1)^2);
    end
end


[rho2d, Teta2d] = calc_cyl_2d(coords2d_r,coords2d_fi);

h2d = zeros(size(coords2d_r));

for i = 1:size(yn,1)
    for j = 1:size(xn,2)
        
        if coords2d_r(i,j) > 50
            h2d(i,j) = NaN;
        elseif sqrt((xn(1,j)+10)^2 + yn(i,1)^2) < 3.95
            h2d(i,j) = NaN;
        else
            h2d(i,j) = A0*R0*(R - rho2d(i,j))/(rho2d(i,j)*(R - R0)) + A1*(R - rho2d(i,j))*cos(Teta2d(i,j))/(R - R0);
        end
        
    end
end


format long

fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Classic_set_gen\randomized_data.txt','r');
formatSpec = ['%f' '%f' '%f'];
sizeMexp = [3 NLearn+NTest];
Mexp = fscanf(fileID,formatSpec,sizeMexp);
fclose(fileID);
Mexp = Mexp';

xexp = Mexp(:,1);
yexp = Mexp(:,2);
hexp = Mexp(:,3);

%solution and experimental data plot
figure('Color','w')
set(gca,'FontSize',12)
a = gradient(0:0.005:0.1);
h = surf(xn, yn, h2d, 'AlphaData',a, 'FaceAlpha',.3);
daspect([1,1,0.03]);
set(h,'edgecolor','r','facecolor',[1 1 1])
xlim([-50 50])
ylim([-50 50])
zlim([0 1.5])
view(154,28)

hold on
plot3(xexp(1:NLearn+NTest),yexp(1:NLearn+NTest),hexp(1:NLearn+NTest),'.k','MarkerSize',20);

xlabel('x(sm)')
ylabel('y(sm)')
zlabel('h(sm)')
legend('h(x,y)', 'experiment',1)

%calculating errors
Arr_calc = zeros(1,NLearn+NTest);
[X, Y] = meshgrid(xn, yn);

for i = 1:NLearn+NTest
    Arr_calc(i) = griddata(X,Y,h2d,xexp(i),yexp(i),'natural');
end

Arr_calc_t = Arr_calc';

err = Arr_calc_t - hexp;
disp(err);
max_err_learn = max(err(1:NLearn));
min_err_learn = min(err(1:NLearn));
if abs(max_err_learn) > abs(min_err_learn)
    max_desc_learn = max_err_learn;
else 
    max_desc_learn = min_err_learn;
end
max_desc_learn
max_err_test = max(err(NLearn+1:NLearn+NTest));
min_err_test = min(err(NLearn+1:NLearn+NTest));
if abs(max_err_test) > abs(min_err_test)
    max_desc_test = max_err_test;
else 
    max_desc_test = min_err_test;
end
max_desc_test
err_sq_learn = sqrt((1/size(Arr_calc_t(1:NLearn,1),1))*sum((Arr_calc_t(1:NLearn,1) - hexp(1:NLearn,1)).*(Arr_calc_t(1:NLearn,1) - hexp(1:NLearn,1))))
err_sq_test = sqrt((1/size(Arr_calc_t(NLearn+1:NLearn+NTest,1),1))*sum((Arr_calc_t(NLearn+1:NLearn+NTest,1) - hexp(NLearn+1:NLearn+NTest,1)).*(Arr_calc_t(NLearn+1:NLearn+NTest,1) - hexp(NLearn+1:NLearn+NTest,1))))

fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Multilayer\res.txt','a');
fprintf(fileID,'%f\t%f\t%f\t%f\t',max_desc_learn,max_desc_test,err_sq_learn,err_sq_test);
fclose(fileID);

%plot discrepancy learning and testing
figure('Color','w')
set(gca,'FontSize',12)
plot3(xexp(1:NLearn),yexp(1:NLearn),err(1:NLearn),'.','MarkerSize',20)
grid on
set(gca,'XTick',-50:25:50)
set(gca,'YTick',-50:25:50)
xlim([-50 50])
ylim([-50 50])
zlim([-0.5 0.5])
set(gca,'ZTick',-0.5:0.1:0.5)
xlabel('x(sm)')
ylabel('y(sm)')
zlabel('d(sm)')
view(154,28)

hold on
plot3(xexp(NLearn+1:NLearn+NTest),yexp(NLearn+1:NLearn+NTest),err(NLearn+1:NLearn+NTest),'r.','MarkerSize',20)

legend('discrepancy learning','discrepancy testing',1)

end

function draw_solution_field_2(A0,A1,A2,R0)

R = 1.0;
Nx = 51;
Ny = 51;
xn = linspace(-50,50,Nx);
yn_t = linspace(-50,50,Ny);
yn = yn_t';
coords2d_r = zeros(size(yn,1),size(xn,2));
coords2d_fi = zeros(size(yn,1),size(xn,2));


for i = 1:size(yn,1)
    for j = 1:size(xn,2)
        coords2d_fi(i,j) = atan2(yn(i,1),xn(1,j));
    end
end
 
for i = 1:size(yn,1)
    for j = 1:size(xn,2)
        coords2d_r(i,j) = sqrt(xn(1,j)^2 + yn(i,1)^2);
    end
end


[rho2d, Teta2d] = calc_cyl_2d(coords2d_r,coords2d_fi);

h2d = zeros(size(coords2d_r));

for i = 1:size(yn,1)
    for j = 1:size(xn,2)
        
        if coords2d_r(i,j) > 50
            h2d(i,j) = NaN;
        elseif sqrt((xn(1,j)+10)^2 + yn(i,1)^2) < 3.95
            h2d(i,j) = NaN;
        else
            h2d(i,j) = A0*R0*(R - rho2d(i,j))/(rho2d(i,j)*(R - R0)) + A1*(R - rho2d(i,j))*cos(Teta2d(i,j))/(R - R0) + A2*R0*(R - rho2d(i,j))*cos(2*Teta2d(i,j))*(1/( (4*R0-3*rho2d(i,j))*(R-R0) ));
        end
        
    end
end


format long

fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Classic_set_gen\randomized_data.txt','r');
formatSpec = ['%f' '%f' '%f'];
sizeMexp = [3 NLearn+NTest];
Mexp = fscanf(fileID,formatSpec,sizeMexp);
fclose(fileID);
Mexp = Mexp';

xexp = Mexp(:,1);
yexp = Mexp(:,2);
hexp = Mexp(:,3);

%solution and experimental data plot
figure('Color','w')
set(gca,'FontSize',12)
a = gradient(0:0.005:0.1);
h = surf(xn, yn, h2d, 'AlphaData',a, 'FaceAlpha',.3);
daspect([1,1,0.03]);
set(h,'edgecolor','r','facecolor',[1 1 1])
xlim([-50 50])
ylim([-50 50])
zlim([0 1.5])
view(154,28)

hold on
plot3(xexp(1:NLearn+NTest),yexp(1:NLearn+NTest),hexp(1:NLearn+NTest),'.k','MarkerSize',20);

xlabel('x(sm)')
ylabel('y(sm)')
zlabel('h(sm)')
legend('h(x,y)', 'experiment',1)

%calculating errors
Arr_calc = zeros(1,NLearn+NTest);
[X, Y] = meshgrid(xn, yn);

for i = 1:NLearn+NTest
    Arr_calc(i) = griddata(X,Y,h2d,xexp(i),yexp(i),'natural');
end

Arr_calc_t = Arr_calc';

err = Arr_calc_t - hexp;
disp(err);
max_err_learn = max(err(1:NLearn));
min_err_learn = min(err(1:NLearn));
if abs(max_err_learn) > abs(min_err_learn)
    max_desc_learn = max_err_learn;
else 
    max_desc_learn = min_err_learn;
end
max_desc_learn
max_err_test = max(err(NLearn+1:NLearn+NTest));
min_err_test = min(err(NLearn+1:NLearn+NTest));
if abs(max_err_test) > abs(min_err_test)
    max_desc_test = max_err_test;
else 
    max_desc_test = min_err_test;
end
max_desc_test
err_sq_learn = sqrt((1/size(Arr_calc_t(1:NLearn,1),1))*sum((Arr_calc_t(1:NLearn,1) - hexp(1:NLearn,1)).*(Arr_calc_t(1:NLearn,1) - hexp(1:NLearn,1))))
err_sq_test = sqrt((1/size(Arr_calc_t(NLearn+1:NLearn+NTest,1),1))*sum((Arr_calc_t(NLearn+1:NLearn+NTest,1) - hexp(NLearn+1:NLearn+NTest,1)).*(Arr_calc_t(NLearn+1:NLearn+NTest,1) - hexp(NLearn+1:NLearn+NTest,1))))

fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Multilayer\res.txt','a');
fprintf(fileID,'%f\t%f\t%f\t%f\t',max_desc_learn,max_desc_test,err_sq_learn,err_sq_test);
fclose(fileID);

%plot discrepancy learning and testing
figure('Color','w')
set(gca,'FontSize',12)
plot3(xexp(1:NLearn),yexp(1:NLearn),err(1:NLearn),'.','MarkerSize',20)
grid on
set(gca,'XTick',-50:25:50)
set(gca,'YTick',-50:25:50)
xlim([-50 50])
ylim([-50 50])
zlim([-0.5 0.5])
set(gca,'ZTick',-0.5:0.1:0.5)
xlabel('x(sm)')
ylabel('y(sm)')
zlabel('d(sm)')
view(154,28)

hold on
plot3(xexp(NLearn+1:NLearn+NTest),yexp(NLearn+1:NLearn+NTest),err(NLearn+1:NLearn+NTest),'r.','MarkerSize',20)

legend('discrepancy learning','discrepancy testing',1)

end

function [N] = NLearn
N = 16;
end

function [N] = NTest
N = 8;
end
