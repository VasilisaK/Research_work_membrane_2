function draw_solution_neuron
clear, clc
 
fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Neuron_set\field.txt','r');
formatSpec = ['%f' '%f' '%f'];
sizeM = [3 Inf];
M = fscanf(fileID,formatSpec,sizeM);
fclose(fileID);
M = M';
     
fileID2 = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Classic_set_gen\randomized_data.txt','r');
sizeMexp = [3 NLearn+NTest];
formatSpec = ['\n%f' '%f' '%f'];
Mexp = fscanf(fileID2,formatSpec,sizeMexp);
fclose(fileID2);
Mexp = Mexp';

%solution data
x = M(:,1);
y = M(:,2);
z = M(:,3);
%experimental data
xexp = Mexp(:,1);
yexp = Mexp(:,2);
hexp = Mexp(:,3);

%grid size
Nx = 51;
Ny = 51;

xi = linspace(min(x),max(x),Nx);
yi = linspace(min(y),max(y),Ny);
[X, Y] = meshgrid(xi, yi);
Z = griddata(x,y,z, X,Y, 'natural');

for i = 1:Nx
    for j = 1:Ny
        if sqrt((X(i,j)+10)^2 + Y(i,j)^2) < 3.95
            Z(i,j) = NaN;
        end
    end
end

%solution and experimental data plot
figure('Color','w')
set(gca,'FontSize',12)
a = gradient(0:0.005:1);
h = surf(X, Y, Z, 'AlphaData', a, 'FaceAlpha', .4);
daspect([1,1,0.03]);
set(h,'edgecolor','r','facecolor',[1 1 1])
xlim([-50 50])
ylim([-50 50])
zlim([-0.1 1.5])
xlabel('x(sm)')
ylabel('y(sm)')
zlabel('h(sm)')
view(154,28)

hold on
plot3(xexp(1:NLearn+NTest),yexp(1:NLearn+NTest),hexp(1:NLearn+NTest),'.k','MarkerSize',20)
legend('h(x,y)', 'experiment', 1)

%calculating errors
Arr_calc = zeros(1,NLearn+NTest);

for i = 1:NLearn+NTest
    Arr_calc(i) = griddata(X,Y,Z,xexp(i),yexp(i),'natural');
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
%max_desc_learn
max_err_test = max(err(NLearn+1:NLearn+NTest));
min_err_test = min(err(NLearn+1:NLearn+NTest));
if abs(max_err_test) > abs(min_err_test)
    max_desc_test = max_err_test;
else 
    max_desc_test = min_err_test;
end
%max_desc_test
err_sq_learn = sqrt((1/size(Arr_calc_t(1:NLearn,1),1))*sum((Arr_calc_t(1:NLearn,1) - hexp(1:NLearn,1)).*(Arr_calc_t(1:NLearn,1) - hexp(1:NLearn,1))));
err_sq_test = sqrt((1/size(Arr_calc_t(NLearn+1:NLearn+NTest,1),1))*sum((Arr_calc_t(NLearn+1:NLearn+NTest,1) - hexp(NLearn+1:NLearn+NTest,1)).*(Arr_calc_t(NLearn+1:NLearn+NTest,1) - hexp(NLearn+1:NLearn+NTest,1))));

%fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Neuron_set\res.txt','w');
fileID = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Neuron_set\res.txt','a');
fprintf(fileID,'%f\t%f\t%f\t%f\t',max_desc_learn,max_desc_test,err_sq_learn,err_sq_test);
fclose(fileID);

%plot discrepancy learning and testing 
figure('Color','w')
set(gca,'FontSize',12)
plot3(xexp(1:NLearn),yexp(1:NLearn),err(1:NLearn),'.','MarkerSize',20)
%daspect([1,1,0.005]);
grid on
set(gca,'XTick',-50:25:50)
set(gca,'YTick',-50:25:50)
set(gca,'ZTick',-0.5:0.1:0.5)
xlim([-50 50])
ylim([-50 50])
zlim([-0.5 0.5])
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