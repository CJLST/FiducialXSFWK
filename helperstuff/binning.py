BINS = {
'mass4l': '|105|160|',
'mass4l_zzfloating': '|105|160|',
#'Nj': '|0|1|2|3|4|14|',
#'pT4l': '|0|10|15|20|25|30|35|40|50|60|75|90|120|150|200|10000|',
'pT4l': '|0|10000|',
'pT4l_kL': '|0|45|80|120|200|1300|',
'rapidity4l': '|0|0.1|0.2|0.3|0.4|0.5|0.6|0.7|0.8|0.9|1|1.15|1.3|1.5|2.5|',
'costhetaZ1': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
'costhetaZ2': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
'phi': '|-3.14159265359|-2.35619449019|-1.57079632679|-0.785398163397|0.0|0.785398163397|1.57079632679|2.35619449019|3.14159265359|',
'phi1': '|-3.14159265359|-2.35619449019|-1.57079632679|-0.785398163397|0.0|0.785398163397|1.57079632679|2.35619449019|3.14159265359|',
'costhetastarZZ': '|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
'massZ1': '|40|60|70|80|84|86|88|90|92|10000|',
'massZ2': '|0|20|22|24|26|27|28|29|30|31|32|34|36|40|45|50|10000|',
#'pTj1': '|-2|30|55|95|200|1300|',
#'pTHj': '|-2|0|30|60|110|1300|',
#'mHj': '|-2|110|180|220|300|400|600|3000|',
#'pTj2': '|-2|30|40|65|90|1300|',
#'mjj': '|-2|0|120|300|1300|',
#'absdetajj': '|-100|0|1.6|3.0|10.0|',
#'dphijj': '|-100|-3.14159265359|-1.5707963267948966|0|1.5707963267948966|3.14159265359|',
#'pTHjj': '|-2|0|20|60|1300|',
'TCjmax': '|0|20|30|40|50|75|10000|',
#'TCjmax': '|-2|0|15|20|30|50|80|1000|',
'TBjmax': '|0|45|60|75|100|200|300|10000|',
#'TBjmax': '|-2|0|30|45|75|150|1000|',
'D0m': '|0.0|0.4|0.5|0.6|0.7|0.8|0.9|1.0|',
'Dcp': '|-0.75|-0.25|-0.1|0|0.1|0.25|0.75|',
'D0hp': '|0.0|0.35|0.4|0.45|0.55|0.65|0.75|1.0|',
'Dint': '|0.0|0.7|0.8|0.9|0.95|1.0|',
'DL1': '|0.0|0.45|0.5|0.6|0.7|1.0|',
'DL1Zg': '|0.0|0.35|0.45|0.5|0.55|0.65|1.0|',
'rapidity4l vs pT4l': '|0|0.5|1.0|2.5| vs |0|40|80|150|1300| / |0|45|120|1300| / |0|45|120|1300|',
#'njets_pt30_eta4p7 vs pT4l': '|0|1|2|20| vs |0|15|30|1300| / |0|60|80|120|1300| / |0|100|170|250|1300|',
'pTj1 vs pTj2': '|-100|30| vs |-100|30| / |30|60| vs |30|60| / |60|350| vs |30|60| / |60|350| vs |60|350|',
'pT4l vs pTHj': '|0|1500| vs |-100|0| / |0|85| vs |0|30| / |85|350| vs |0|45| / |0|85| vs |30|350| / |85|350| vs |45|350|',
'massZ1 vs massZ2': '|40|85| vs |12|35| / |40|70| vs |35|65| / |70|120| vs |35|65| / |85|120| vs |30|35| / |85|120| vs |24|30| / |85|120| vs |12|24|',
'TCjmax vs pT4l': '|-2|15|25|40|1300| vs |0|15|30|45|70|120|1300| / |0|120|1300| / |0|120|1300| / |0|200|1300|',


'pTj1': '|-100|30|45|60|90|150|10000|',
'pTj2': '|-100|30|40|50|10000|',
'Nj': '|0|1|2|3|4|10|',
'mjj': '|-100|0|60|90|120|180|240|400|800|10000|',
'absdetajj': '|-100|0|0.5|1|1.5|2|3|4|5|100|',
'dphijj': '|-100|-3.14159265359|-1.5707963267948966|0|1.5707963267948966|3.14159265359|',
'pTHj': '|-100|0|10|20|30|40|50|60|75|100|10000|',
'mHj': '|-100|0|175|200|225|250|275|300|350|400|500|600|800|10000|',
'pTHjj': '|-100|0|15|30|45|60|10000|',

    
'absdetajj vs mjj': '|-100|0| vs |-100|0| / |0|3| vs |0|450| / |0|3| vs |450|10000| / |3|10| vs |0|450| / |3|10| vs |450|10000|',
'pT4l vs Nj': '|0|30| vs |0|1| / |30|10000| vs |0|1| / |0|60| vs |1|2| / |60|10000| vs |1|2| / |0|120| vs |2|10| / |120|10000| vs |2|10|',
'rapidity4l vs pT4l': '|0|1| vs |0|30| / |0|1| vs |30|10000| / |1|2| vs |0|60| / |1|2| vs |60|10000| / |2|10| vs |0|120| / |2|10| vs |120|10000|',
'pT4l vs pTj1': '|0|30| vs |0|45| / |0|30| vs |45|10000| / |30|60| vs |0|90| / |30|60| vs |90|10000| / |60|10000| vs |0|120| / |60|10000| vs |120|10000|',
    
}

def binning(var):
    obsBins_input = BINS[var]
    if not 'vs' in obsBins_input: #It is not a double-differential analysis
        obs_bins = {0:(obsBins_input.split("|")[1:(len(obsBins_input.split("|"))-1)]),1:['0','inf']}[obsBins_input=='inclusive']
        obs_bins = [float(i) for i in obs_bins] #Convert a list of str to a list of float
        doubleDiff = False
        print ('It is a single-differential measurement, binning', obs_bins)
    else: #It is a double-differential analysis
        doubleDiff = True
        # The structure of obs_bins is:
        # index of the dictionary is the number of the bin
        # [obs_bins_low, obs_bins_high, obs_bins_low_2nd, obs_bins_high_2nd]
        # The first two entries are the lower and upper bound of the first variable
        # The second two entries are the lower and upper bound of the second variable
        if obsBins_input.count('vs')==1 and obsBins_input.count('/')>=1: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|'
            obs_bins_tmp = obsBins_input.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250| / |0|10|20|80|250| / |0|20|90|250| / |0|25|250|']
            obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
            obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
            obs_bins_tmp = obs_bins_tmp[1].split(' / ') #['|0|10|20|45|90|250|', '|0|10|20|80|250|', '|0|20|90|250|', '|0|25|250|']
            obs_bins_2nd = {}
            for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_2nd {0: ['0', '10', '20', '45', '90', '250'], 1: ['0', '10', '20', '80', '250'], 2: ['0', '20', '90', '250'], 3: ['0', '25', '250']}
                obs_bins_2nd[i] = obs_bins_tmp[i].split('|')[1:len(obs_bins_tmp[i].split('|'))-1]
                obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
            obs_bins = {}
            k = 0 #Bin index
            for i in range(len(obs_bins_1st)-1):
                for j in range(len(obs_bins_2nd[i])-1):
                    obs_bins[k] = []
                    obs_bins[k].append(obs_bins_1st[i])
                    obs_bins[k].append(obs_bins_1st[i+1])
                    obs_bins[k].append(obs_bins_2nd[i][j])
                    obs_bins[k].append(obs_bins_2nd[i][j+1])
                    k +=1
        elif obsBins_input.count('vs')>1 and obsBins_input.count('/')>1: #Situation like this one '|50|80| vs |10|30| / |50|80| vs |30|60| / |80|110| vs |10|25| / |80|110| vs |25|30|'
            obs_bins_tmp = obsBins_input.split(' / ') #['|50|80| vs |10|30|', '|50|80| vs |30|60|', '|80|110| vs |10|25|', '|80|110| vs |25|30|']
            obs_bins_1st={}
            obs_bins_2nd={}
            obs_bins={}
            for i in range(len(obs_bins_tmp)): #At the end of the loop -> obs_bins_1st {0: ['50', '80'], 1: ['50', '80'], 2: ['80', '110'], 3: ['80', '110']} and obs_bins_2nd {0: ['10', '30'], 1: ['30', '60'], 2: ['10', '25'], 3: ['25', '30']}
                obs_bins_tmp_bis = obs_bins_tmp[i].split(' vs ')
                obs_bins_1st[i] = obs_bins_tmp_bis[0].split('|')[1:len(obs_bins_tmp_bis[0].split('|'))-1]
                obs_bins_1st[i] = [float(j) for j in obs_bins_1st[i]] #Convert a list of str to a list of float
                obs_bins_2nd[i] = obs_bins_tmp_bis[1].split('|')[1:len(obs_bins_tmp_bis[1].split('|'))-1]
                obs_bins_2nd[i] = [float(j) for j in obs_bins_2nd[i]] #Convert a list of str to a list of float
                obs_bins[i] = []
                obs_bins[i].append(obs_bins_1st[i][0])
                obs_bins[i].append(obs_bins_1st[i][1])
                obs_bins[i].append(obs_bins_2nd[i][0])
                obs_bins[i].append(obs_bins_2nd[i][1])
        elif obsBins_input.count('vs')==1 and obsBins_input.count('/')==0: #Situation like this one '|0|1|2|3|20| vs |0|10|20|45|90|250|'
            obs_bins_tmp = obsBins_input.split(" vs ") #['|0|1|2|3|20|', '|0|10|20|45|90|250|']
            obs_bins_1st = obs_bins_tmp[0].split('|')[1:len(obs_bins_tmp[0].split('|'))-1] #['0', '1', '2', '3', '20']
            obs_bins_1st = [float(i) for i in obs_bins_1st] #Convert a list of str to a list of float
            obs_bins_2nd = obs_bins_tmp[1].split('|')[1:len(obs_bins_tmp[1].split('|'))-1] #['0', '10', '20', '45', '90', '250']
            obs_bins_2nd = [float(i) for i in obs_bins_2nd] #Convert a list of str to a list of float
            obs_bins = {}
            k = 0 #Bin index
            for i in range(len(obs_bins_1st)-1):
                for j in range(len(obs_bins_2nd)-1):
                    obs_bins[k] = []
                    obs_bins[k].append(obs_bins_1st[i])
                    obs_bins[k].append(obs_bins_1st[i+1])
                    obs_bins[k].append(obs_bins_2nd[j])
                    obs_bins[k].append(obs_bins_2nd[j+1])
                    k +=1
        else:
            print ('Problem in the definition of the binning')
            quit()
        print ('It is a double-differential measurement, binning for the 1st variable', obs_bins_1st, 'and for the 2nd variable', obs_bins_2nd)
        print (obs_bins)
    return obs_bins, doubleDiff

def binning_v2(var):
    obs_bins, doubleDiff = binning(var)
    return obs_bins
