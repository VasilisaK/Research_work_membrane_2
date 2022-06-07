Xg = [];
Yg = [];
G = [];

with open('C:/Users/VasilisaK/Desktop/Current/NIR_2/Classic_set_gen/randomized_data.txt') as f:
    
    for line in f:
        s = line.strip().split()
        Xg.append(float(s[0]))
        Yg.append(float(s[1]))
        G.append(float(s[2]))
        
with open('wolfram_data.txt', 'w') as f:
    
    counter = 0
    clines = 0;
    for x in Xg:
        counter += 1
        s = 'Xg[' + str(counter) + '] = ' + str(x) + ';\n'
        f.write(s)
        clines += 1
        if clines == 16:
            break
        
    f.write('\n')
    
    counter = 0
    clines = 0;
    for y in Yg:
        counter += 1
        s = 'Yg[' + str(counter) + '] = ' + str(y) + ';\n'
        f.write(s)
        clines += 1
        if clines == 16:
            break
         
    f.write('\n')
    
    counter = 0
    clines = 0;
    for g in G:
        counter += 1
        s = 'G[' + str(counter) + '] = ' + str(g) + ';\n'
        f.write(s)
        clines += 1
        if clines == 16:
            break
