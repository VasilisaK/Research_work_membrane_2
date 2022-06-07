function generate_sample
clear, clc

numb = NLearn;
numbers_learn = randsample(1:NLearn+NTest,numb);
numbers_extra = zeros(NLearn+NTest,1);
numbers_final = zeros(NLearn+NTest,1);

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

for i = 1:numb
    numbers_extra(numbers_learn(i)) = 10e8;
end

for i = 1:numb
    numbers_final(i) = numbers_learn(i);
end

j = numb;
for i = 1:NLearn+NTest
    if numbers_extra(i) < 10e8
        j = j + 1;
        numbers_final(j) = i;        
    end
end

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

fileID3 = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Classic_set_gen\numbers.txt','w');
fprintf(fileID3,'%i\t',numbers_final);
fprintf(fileID3,'\n');
fclose(fileID3);

fileID4 = fopen('C:\Users\VasilisaK\Desktop\Current\NIR_2\Classic_set_gen\numbers_combinations.txt','a');
fprintf(fileID4,'%i\t',numbers_final);
fprintf(fileID4,'\n');
fclose(fileID4);

end

function [N] = NLearn
N = 16;
end

function [N] = NTest
N = 8;
end