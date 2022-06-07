function make_sample
clear, clc

% reading experimental data
fileID = fopen('experimental_data.txt','r');
formatSpec = ['%f' '%f' '%f'];
sizeM = [3 NLearn+NTest];
M = fscanf(fileID,formatSpec,sizeM);
fclose(fileID);
M = M';

x_raw = M(1:NLearn+NTest,1);
y_raw = M(1:NLearn+NTest,2);
hexp_raw = M(1:NLearn+NTest,3);

x = zeros(NLearn+NTest,1);
y = zeros(NLearn+NTest,1);
hexp = zeros(NLearn+NTest,1);

fileID2 = fopen('numbers.txt','r');
numbers_final = fscanf(fileID2,'%f');
fclose(fileID2);

for i = 1:NLearn+NTest
    x(i) = x_raw(numbers_final(i));
    y(i) = y_raw(numbers_final(i));
    hexp(i) = hexp_raw(numbers_final(i));
end

% writing randomized experimental data
fileID2 = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Classic_set_gen\randomized_data.txt','w');
for i = 1:NLearn+NTest
    fprintf(fileID,'%.15f %.15f %.15f\n',x(i),y(i),hexp(i));
end
fclose(fileID2);

end

function [N] = NLearn
N = 16;
end

function [N] = NTest
N = 8;
end