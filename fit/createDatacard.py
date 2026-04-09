import os,sys
from numpy import array, float32 # spencer

def fixJes(jesnp, jes_evts_noWeight):
    if jes_evts_noWeight <= 50:
        return '- '
     # Cases where jes are '-' or single value
    if (len(jesnp)==1) | ('/' not in jesnp):
        return jesnp+' '
    else:
        jesnp_tmp = jesnp.split('/')
        jesnp_tmp_dn = jesnp_tmp[0]
        jesnp_tmp_up = jesnp_tmp[1]
        if ((jesnp_tmp_dn=='1.0') & (jesnp_tmp_up!='1.0')):
            return jesnp_tmp_up+' '
        elif ((jesnp_tmp_dn!='1.0') & (jesnp_tmp_up=='1.0')):
            return jesnp_tmp_dn+'/'+ str(abs(round(2-float(jesnp_tmp_dn),3)))+' '
        elif ((float(jesnp_tmp_dn) > 1) & (float(jesnp_tmp_up) > 1)):
            return jesnp_tmp_up + '/' + str(abs(round(2-float(jesnp_tmp_up),3)))+' '
        elif ((float(jesnp_tmp_dn) < 1) & (float(jesnp_tmp_up) < 1)):
            return jesnp_tmp_dn + '/' + str(abs(round(2-float(jesnp_tmp_dn),3)))+' '
        elif (float(jesnp_tmp_dn) > float(jesnp_tmp_up)):
            '''
            if 1.X/0.X cases
            return a single NP (symmetric) corresponding to largest variation
            return the variation always as 1.X
            '''
            if((float(jesnp_tmp_dn)-1) > (1-float(jesnp_tmp_up))):
                return jesnp_tmp_dn+' '
            else:
                return str(1+(1-float(jesnp_tmp_up)))+' '
        else:
            '''
            if dn/up: correct jes, use it
            '''
            return jesnp+' '

def createDatacard(obsName, channel, nBins, obsBin, observableBins, physicalModel, year, nData, jes, lowerBound, upperBound, yearSetting):
    # Name of the bin (aFINALSTATE_ recobinX)
    if(channel == '4mu'): channelNumber = 1
    if(channel == '4e'): channelNumber = 2
    if(channel == '2e2mu'): channelNumber = 3

    # ZZfloating
    if 'zzfloating' in obsName: zzfloating = True
    else: zzfloating = False

    if '_' in obsName and not 'floating' in obsName and not 'kL' in obsName and not obsName == 'Nj': #it means it is a double differential measurement
        _recobin = str(observableBins[obsBin][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin][3]).replace('.', 'p').replace('-','m')
    else:
        _recobin = str(observableBins[obsBin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin+1]).replace('.', 'p').replace('-','m')
        if int(observableBins[obsBin+1]) > 1000:
            _recobin = 'GT'+str(int(observableBins[obsBin]))

    if physicalModel == 'v3':
        _obsName = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'pTj1', 'Nj': 'Nj'}
        if zzfloating:
            if obsName.rsplit('_', 1)[0] not in _obsName: 
                _obsName[obsName] = obsName
            else:
                _obsName[obsName] = _obsName[obsName.rsplit('_', 1)[0]]+'_zzfloating'
        else:
            if obsName not in _obsName:
                _obsName[obsName] = obsName 
            else:
                _obsName[obsName] = _obsName[obsName]
        binName = 'hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel
        processName = 'smH_' + _obsName[obsName]
    else:
        _obsName = {}
        _obsName[obsName] = obsName
        binName = 'a'+str(channelNumber)+'_rerfcobin'+str(obsBin)
        processName = 'trueH'+channel+'Bin'

    # Background expectations
    sys.path.append('../inputs')
    #_temp = __import__('inputs_bkgTemplate_'+obsName, globals(), locals(), ['expected_yield'], -1)
    _temp = __import__('inputs_bkgTemplate_'+obsName, globals(), locals(), ['expected_yield'], 0) # spencer 
    expected_yield = _temp.expected_yield
    #_temp = __import__('inputs_bkg_'+obsName+'_'+year, globals(), locals(), ['fractionsBackground'], -1)
    _temp = __import__('inputs_bkg_'+obsName+'_'+year, globals(), locals(), ['fractionsBackground'], 0) # spencer 
    fractionsBackground = _temp.fractionsBackground
    if jes:

        sys.path.append('../coefficients/JES')
        #_temp = __import__('JESNP_'+obsName, globals(), locals(), ['JESNP'], -1)
        if 'zzfloating' in obsName: _temp = __import__('JESNP_'+obsName.rsplit('_', 1)[0]+'_'+year, globals(), locals(), ['JESNP'], 0) # spencer
        else: _temp = __import__('JESNP_'+obsName+'_'+year, globals(), locals(), ['JESNP'], 0) # spencer
        jesnp = _temp.JESNP
        #_temp = __import__('JESNP_evts_'+obsName, globals(), locals(), ['evts_noWeight'], -1)
        if 'zzfloating' in obsName: _temp = __import__('JESNP_evts_'+obsName.rsplit('_', 1)[0]+'_'+year, globals(), locals(), ['evts_noWeight'], 0) # spencer
        else: _temp = __import__('JESNP_evts_'+obsName+'_'+year, globals(), locals(), ['evts_noWeight'], 0) # spencer
        jes_evts_noWeight = _temp.evts_noWeight

        if year == "2023preBPix":
            year_jes = "2023"
        elif year == "2023postBPix":
            year_jes = "2023BPix"
        else:
            year_jes = year 

        jesNames = ['Absolute', 'Absolute_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavorQCD', 'HF', 'HF_year', 'RelativeBal', 'RelativeSample_year']
        jesNames_datacard = [j.replace('year',year_jes) for j in jesNames] # The name of the nuisance in the datacard should have the correspoding year
        sys.path.remove('../coefficients/JES')
    sys.path.remove('../inputs')

    # lumi
    if yearSetting == 'Full':
        lumi = {}
        lumi['2016'] = '1.01'
        lumi['2017'] = '1.02'
        lumi['2018'] = '1.015'
        lumi_corr_16_17_18 = {}
        lumi_corr_16_17_18['2016'] = '1.006'
        lumi_corr_16_17_18['2017'] = '1.009'
        lumi_corr_16_17_18['2018'] = '1.02'
        lumi_corr_17_18 = {}
        lumi_corr_17_18['2017'] = '1.006'
        lumi_corr_17_18['2018'] = '1.002'
    elif yearSetting == 'Run3':
        # RUN III correlated luminosity scheme (2022/23/24)
        lumi = {
        '2022':'1.0138',
        '2022EE':'1.0138',
        '2023preBPix':'1.0017',
        '2023postBPix':'1.0017',
        '2024':'1.0020',
        }
        # Additional reduced nuisances for Run-3:
        # lumi_13p6TeV_2324 applies to 2023 & 2024
        # lumi_13p6TeV_2024 applies to 2024 only
        lumi_extra = {
        'lumi_13p6TeV_2324': {'2023preBPix': '1.0127', '2023postBPix': '1.0127', '2024': '1.0068'},
        'lumi_13p6TeV_2024': {'2024': '1.0144'}
        }
    else:
        lumi = {}
        lumi['2016'] = '1.026'
        lumi['2017'] = '1.025'
        lumi['2018'] = '1.023'

        # RUN III: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiRecommendationsRun3
        # 2024: https://indico.cern.ch/event/1608377/contributions/6777221/attachments/3166738/5628831/2024%20pp%20plots%20and%20detector%20information%20(2).pdf 
        lumi['2022'] = '1.014'
        lumi['2022EE'] = '1.014'
        lumi['2023preBPix'] = '1.013'
        lumi['2023postBPix'] = '1.013'
        lumi['2024'] = '1.0161' 


    # Lepton efficiency
    # Values taken from:
    # https://indico.cern.ch/event/1125278/contributions/4723319/attachments/2420259/4142589/LepSyst.pdf
    # Computed with:
    # https://github.com/asculac/LeptonSystematics
    eff_mu = {}
    eff_mu['2016_2e2mu'] = '0.986/1.007'
    eff_mu['2016_4mu'] = '0.981/1.01'
    eff_mu['2017_2e2mu'] = '0.986/1.006'
    eff_mu['2017_4mu'] = '0.981/1.009'
    eff_mu['2018_2e2mu'] = '0.986/1.006'
    eff_mu['2018_4mu'] = '0.981/1.008'

    # Latest updated values for Run3 (Muons efficiencies, Trigger efficiency is set to 0 here)
    eff_mu['2022_2e2mu'] = '0.989/1.011' #ok
    eff_mu['2022_4mu'] = '0.986/1.014' #ok

    eff_mu['2022EE_2e2mu'] = '0.993/1.007' #ok
    eff_mu['2022EE_4mu'] = '0.992/1.008'#ok

    eff_mu['2023preBPix_2e2mu'] = '0.995/1.005' # ok
    eff_mu['2023preBPix_4mu'] = '0.994/1.006' # ok

    eff_mu['2023postBPix_2e2mu'] = '0.994/1.006' # ok
    eff_mu['2023postBPix_4mu'] = '0.992/1.008' # ok

    eff_mu['2024_2e2mu'] = '0.997/1.003' # ok
    eff_mu['2024_4mu'] = '0.997/1.003' # ok
    
    eff_e = {}
    eff_e['2016_2e2mu'] = '0.934/1.062'
    eff_e['2016_4e'] = '0.891/1.093'
    eff_e['2017_2e2mu'] = '0.953/1.043'
    eff_e['2017_4e'] = '0.915/1.064'
    eff_e['2018_2e2mu'] = '0.95/1.052'
    eff_e['2018_4e'] = '0.905/1.077'

    #CORRELATED

    eff_e['2022_4e'] = '0.884/1.103' #perfect agreement w/ our measurement
    eff_e['2022_2e2mu'] = '0.928/1.069' #perfect agreement w/ our measurement
    eff_e['2022EE_4e'] =  '0.897/1.088' #perfect agreement w/ our measurement
    eff_e['2022EE_2e2mu'] = '0.938/1.059' #perfect agreement w/ our measurement

    eff_e['2023preBPix_2e2mu'] = '0.882/1.116' # ok
    eff_e['2023preBPix_4e'] = '0.814/1.177' # ok
    eff_e['2023postBPix_2e2mu'] = '0.852/1.146' # ok
    eff_e['2023postBPix_4e'] = '0.778/1.215' # ok

    eff_e['2024_2e2mu'] = '0.933/1.064'# ok
    eff_e['2024_4e'] = '0.887/1.098'# ok

    #DECORRELATED
    eff_e_reco_stat = {}
    eff_e_reco_syst = {}
    eff_e_id_stat = {}
    eff_e_id_syst = {}

    # 2022
    eff_e_reco_stat['2022_4e'] = '0.947/1.053' 
    eff_e_reco_stat['2022_2e2mu'] = '0.964/1.036'

    eff_e_reco_syst['2022_4e'] = '0.953/1.047' 
    eff_e_reco_syst['2022_2e2mu'] = '0.967/1.033'

    eff_e_id_stat['2022_4e'] = '0.957/1.043' 
    eff_e_id_stat['2022_2e2mu'] = '0.971/1.029' 

    eff_e_id_syst['2022_4e'] = '0.986/1.014' 
    eff_e_id_syst['2022_2e2mu'] = '0.991/1.009'

    # 2022EE
    eff_e_reco_stat['2022EE_4e'] = '0.946/1.054' 
    eff_e_reco_stat['2022EE_2e2mu'] = '0.964/1.036'

    eff_e_reco_syst['2022EE_4e'] = '0.9503/1.0497' 
    eff_e_reco_syst['2022EE_2e2mu'] = '0.966/1.034'

    eff_e_id_stat['2022EE_4e'] = '0.978/1.022' 
    eff_e_id_stat['2022EE_2e2mu'] = '0.985/1.015' 

    eff_e_id_syst['2022EE_4e'] = '0.982/1.018' 
    eff_e_id_syst['2022EE_2e2mu'] = '0.988/1.012'

    # 2023preBPix
    eff_e_reco_stat['2023preBPix_4e'] = '0.945/1.055'
    eff_e_reco_stat['2023preBPix_2e2mu'] = '0.964/1.036'

    eff_e_reco_syst['2023preBPix_4e'] = '0.841/1.159' 
    eff_e_reco_syst['2023preBPix_2e2mu'] = '0.897/1.103'

    eff_e_id_stat['2023preBPix_4e'] = '0.962/1.038' 
    eff_e_id_stat['2023preBPix_2e2mu'] = '0.976/1.024'

    eff_e_id_syst['2023preBPix_4e'] = '0.990/1.010' 
    eff_e_id_syst['2023preBPix_2e2mu'] = '0.994/1.006'

    # 2023postBPix
    eff_e_reco_stat['2023postBPix_4e'] = '0.924/1.076' 
    eff_e_reco_stat['2023postBPix_2e2mu'] = '0.949/1.051'

    eff_e_reco_syst['2023postBPix_4e'] = '0.838/1.162'
    eff_e_reco_syst['2023postBPix_2e2mu'] = '0.890/1.110'

    eff_e_id_stat['2023postBPix_4e'] = '0.941/1.059' 
    eff_e_id_stat['2023postBPix_2e2mu'] = '0.958/1.042' 

    eff_e_id_syst['2023postBPix_4e'] = '0.972/1.028' 
    eff_e_id_syst['2023postBPix_2e2mu'] = '0.981/1.019'

    # 2024 (to cross check numbers)
    eff_e_reco_stat['2024_4e'] = '0.944/1.056' 
    eff_e_reco_stat['2024_2e2mu'] = '0.962/1.038'

    eff_e_reco_syst['2024_4e'] = '0.949/1.051' 
    eff_e_reco_syst['2024_2e2mu'] = '0.965/1.035'

    eff_e_id_stat['2024_4e'] = '0.978/1.022' 
    eff_e_id_stat['2024_2e2mu'] = '0.983/1.017' 

    eff_e_id_syst['2024_4e'] = '0.982/1.018' 
    eff_e_id_syst['2024_2e2mu'] = '0.987/1.013'


    #Trigger (from Run 2)
    trig_mu = {}

    trig_mu['2022_4mu'] = '0.983/1.001'
    trig_mu['2022_2e2mu'] = '0.988/1.002'

    trig_mu['2022EE_4mu'] =  '0.983/1.001'
    trig_mu['2022EE_2e2mu'] = '0.988/1.002'

    trig_mu['2023preBPix_2e2mu'] = '0.988/1.002'
    trig_mu['2023preBPix_4mu'] = '0.983/1.001'
    trig_mu['2023postBPix_2e2mu'] = '0.988/1.002'
    trig_mu['2023postBPix_4mu'] = '0.983/1.001'

    trig_mu['2024_2e2mu'] = '0.988/1.002'
    trig_mu['2024_4mu'] = '0.983/1.001'

    trig_e = {}
    trig_e['2022_4e'] = '0.945/1.011' 
    trig_e['2022_2e2mu'] = '0.979/1.002'

    trig_e['2022EE_4e'] =  '0.945/1.011'
    trig_e['2022EE_2e2mu'] = '0.979/1.002'

    trig_e['2023preBPix_2e2mu'] = '0.979/1.002'
    trig_e['2023preBPix_4e'] = '0.945/1.011'
    trig_e['2023postBPix_2e2mu'] = '0.979/1.002'
    trig_e['2023postBPix_4e'] = '0.945/1.011'

    trig_e['2024_2e2mu'] = '0.979/1.002'
    trig_e['2024_4e'] = '0.945/1.011'


    # ZX
    ZX = {}
    ZX['2016_2e2mu'] = '0.756295/1.25114'
    ZX['2016_4e'] = '0.598182/1.43059'
    ZX['2016_4mu'] = '0.694678/1.30555'
    ZX['2017_2e2mu'] = '0.765736/1.236385'
    ZX['2017_4e'] = '0.646521/1.36398'
    ZX['2017_4mu'] = '0.694063/1.30623'
    ZX['2018_2e2mu'] = '0.7660052/1.235647'
    ZX['2018_4e'] = '0.650486/1.35893'
    ZX['2018_4mu'] = '0.69554/1.30465'

    # 2022 values taken from https://indico.cern.ch/event/1360904/contributions/5896503/attachments/2832114/4948339/HZZmeeting_050424.pdf (slide 14)
    ZX['2022_2e2mu'] = '0.724/1.263'
    ZX['2022_4e'] = '0.495/1.451'
    ZX['2022_4mu'] = '0.677/1.321'

    ZX['2022EE_2e2mu'] = '0.748/1.245'
    ZX['2022EE_4e'] = '0.575/1.398'
    ZX['2022EE_4mu'] = '0.690/1.310'

    ZX['2023preBPix_2e2mu'] = '0.734/1.256' # spencer
    ZX['2023preBPix_4e'] = '0.516/1.438' # spencer
    ZX['2023preBPix_4mu'] = '0.685/1.314' # spencer

    ZX['2023postBPix_2e2mu'] = '0.709/1.272' # spencer
    ZX['2023postBPix_4e'] = '0.467/1.470' # spencer
    ZX['2023postBPix_4mu'] = '0.667/1.329' # spencer 

    ZX['2024_2e2mu'] = '0.766/1.232' # spencer
    ZX['2024_4e'] = '0.654/1.339' # spencer
    ZX['2024_4mu'] = '0.696/1.303' # spencer
    
    
    # -------------------------------------------------------------------------------------------------


    print('OPENING ../datacard/datacard_'+year+'/hzz4l_'+channel+'S_13TeV_xs_'+_obsName[obsName]+'_bin'+str(obsBin)+'_'+physicalModel+'.txt')
    file = open('../datacard/datacard_'+year+'/hzz4l_'+channel+'S_13TeV_xs_'+_obsName[obsName]+'_bin'+str(obsBin)+'_'+physicalModel+'.txt', 'w+')

    file.write('imax 1 \n')
    file.write('jmax * \n')
    file.write('kmax * \n')

    file.write('------------ \n')

    file.write('shapes * * hzz4l_'+channel+'S_13TeV_xs_SM_125_'+_obsName[obsName]+'_'+physicalModel+'.Databin'+str(obsBin)+'.root w:$PROCESS\n')

    file.write('------------ \n')

    file.write('bin '+binName+'\n')
    file.write('observation '+str(nData)+'\n')

    file.write('------------ \n')
    file.write('## mass window ['+str(lowerBound)+','+str(upperBound)+']\n')
    file.write('bin ')
    # for i in range(nBins+5): # In addition to the observableBins, there are OutsideAcceptance, fakeH, bkg_ggzz, bkg_qqzz, bkg_zjets
    for i in range(nBins+5):
        file.write(binName+' ')
    file.write('\n')
    file.write('process ')
    if physicalModel == 'v3':
        for i in range(nBins):
            if '_' in obsName and not 'floating' in obsName and not 'kL' in obsName and not obsName == 'Nj':
                file.write(processName+'_'+str(observableBins[i][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][3]).replace('.', 'p').replace('-','m')+' ')
            elif observableBins[i+1] > 1000:
                file.write(processName+'_GT'+str(int(observableBins[i]))+' ')
            else:
                file.write(processName+'_'+str(observableBins[i]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i+1]).replace('.', 'p').replace('-','m')+' ')
        file.write('OutsideAcceptance nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    else:
        for i in range(nBins):
            file.write(processName+str(i)+' ')
        file.write('out_trueH fakeH bkg_qqzz bkg_ggzz bkg_zjets')
    #file.write('nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    file.write('\n')
    file.write('process ')
    for i in range(nBins):
        file.write('-'+str(i+1)+' ')
    # file.write('1 2 3 4 5')
    file.write('1 2 3 4 5')
    file.write('\n')
    file.write('rate ')
    # for i in range(nBins+2): # In addition to the observableBins, there are OutsideAcceptance, fakeH
    for i in range(nBins+2):
        file.write('1.0 ')
    if zzfloating:
        file.write('1 1 '+str(expected_yield[year,'ZX',channel])+'\n')
    else:
        print(year)
        print(expected_yield[year,'qqzz',channel])
        file.write(str(expected_yield[year,'qqzz',channel])+' '+str(expected_yield[year,'ggzz',channel])+' '+str(expected_yield[year,'ZX',channel])+'\n')
    file.write('------------ \n')

    if zzfloating:
        # rateParam qqZZ floating
        if physicalModel == 'v3':
            #min_range = 0
            #file.write('zz_norm_'+str(obsBin)+' rateParam '+binName+' bkg_*zz '+str(expected_yield['ZZ_'+str(obsBin)])+' ['+str(min_range)+','+str(expected_yield['ZZ_'+str(obsBin)]+100)+']\n')
            file.write('zz_norm_'+str(obsBin)+' rateParam '+binName+' bkg_*zz '+str(expected_yield['ZZ_'+str(obsBin)])+' ['+str(expected_yield['ZZ_'+str(obsBin)]*0.5)+','+str(expected_yield['ZZ_'+str(obsBin)]*1.5)+']\n')

        elif physicalModel == 'v2':
            if channel == '2e2mu':
                min_range = 0
            else:
                min_range = expected_yield['ZZ_'+channel]-100
            file.write('zz_norm_'+str(obsBin)+'_'+channel+' rateParam '+binName+' bkg_*zz '+str(expected_yield['ZZ_'+channel])+' ['+str(min_range)+','+str(expected_yield['ZZ_'+channel]+100)+']\n')

    if yearSetting == 'Full':
        if zzfloating:
            # lumi_uncorrelated
            file.write('lumi_13TeV_'+year+' lnN ')
            for i in range(nBins+2): # signals + out + fake
                file.write(lumi[year]+' ')
            file.write('- - -\n') # qqzz + ggzz + ZX
            # lumi_correlated_16_17_18
            file.write('lumi_13TeV_correlated lnN ')
            for i in range(nBins+2): # signals + out + fake
                file.write(lumi_corr_16_17_18[year]+' ')
            file.write('- - -\n') # qqzz + ggzz + ZX
            # lumi_correlated_17_18
            if year == '2017' or year == '2018':
                file.write('lumi_13TeV_1718 lnN ')
                for i in range(nBins+2): # signals + out + fake
                    file.write(lumi_corr_17_18[year]+' ')
                file.write('- - -\n') # qqzz + ggzz + ZX
        else:
            # lumi_uncorrelated
            file.write('lumi_13TeV_'+year+' lnN ')
            for i in range(nBins+4): # All except ZX
                file.write(lumi[year]+' ')
            file.write('-\n') # ZX
            # lumi_correlated_16_17_18
            file.write('lumi_13TeV_correlated lnN ')
            for i in range(nBins+4): # All except ZX
                file.write(lumi_corr_16_17_18[year]+' ')
            file.write('-\n') # ZX
            # lumi_correlated_17_18
            if year == '2017' or year == '2018':
                file.write('lumi_13TeV_1718 lnN ')
                for i in range(nBins+4): # All except ZX
                    file.write(lumi_corr_17_18[year]+' ')
                file.write('-\n') # ZX

    elif yearSetting == 'Run3':
        if zzfloating:
            # correlated Run-3 lumi nuisances
            # lumi_13p6TeV_222324 (2022,2023,2024)
            file.write('lumi_13p6TeV_222324 lnN ')
            for i in range(nBins+2): 
                # All except ZX
                file.write(lumi[year]+' ')
            file.write('- - -\n')

            # lumi_13p6TeV_2324 if applicable
            if year in lumi_extra['lumi_13p6TeV_2324']:
                file.write('lumi_13p6TeV_2324 lnN ')
                for i in range(nBins+2):
                    file.write(lumi_extra['lumi_13p6TeV_2324'][year]+' ')
                file.write('- - -\n')

            # lumi_13p6TeV_2024 if applicable
            if year in lumi_extra['lumi_13p6TeV_2024']:
                file.write('lumi_13p6TeV_2024 lnN ')
                for i in range(nBins+2):
                    file.write(lumi_extra['lumi_13p6TeV_2024'][year]+' ')
                file.write('- - -\n')
        else:
            # correlated Run-3 lumi nuisances
            # lumi_13p6TeV_222324 (2022,2023,2024)
            file.write('lumi_13p6TeV_222324 lnN ')
            for i in range(nBins+4): 
                # All except ZX
                file.write(lumi[year]+' ')
            file.write('-\n')

            # lumi_13p6TeV_2324 if applicable
            if year in lumi_extra['lumi_13p6TeV_2324']:
                file.write('lumi_13p6TeV_2324 lnN ')
                for i in range(nBins+4):
                    file.write(lumi_extra['lumi_13p6TeV_2324'][year]+' ')
                file.write('-\n')

            # lumi_13p6TeV_2024 if applicable
            if year in lumi_extra['lumi_13p6TeV_2024']:
                file.write('lumi_13p6TeV_2024 lnN ')
                for i in range(nBins+4):
                    file.write(lumi_extra['lumi_13p6TeV_2024'][year]+' ')
                file.write('-\n')

    else:

        # lumi
        if year == '2022' or year == '2022EE':
            file.write('lumi_13p6TeV_2022 lnN ')

        elif year == '2023preBPix' or year == '2023postBPix':
            file.write('lumi_13p6TeV_2023 lnN ')

        elif year == '2024':
            file.write('lumi_13p6TeV_2024 lnN ')
        else:
            file.write('lumi_13TeV_'+year+' lnN ')

        if zzfloating:
            for i in range(nBins+2): # All except ZX
                file.write(lumi[year]+' ')
            file.write('- - -\n') # ZX
        else:
            for i in range(nBins+4): # All except ZX
                file.write(lumi[year]+' ')
            file.write('-\n') # ZX

            
    # Lepton efficiency

    if channel == '4mu' or channel == '2e2mu':
        file.write('CMS_eff_m lnN ')
        for i in range(nBins+4): # All except ZX
            file.write(eff_mu[year+'_'+channel]+' ')
        file.write('-\n') 

        if year == "2023preBPix":
            file.write('CMS_eff_m_trigger_2023 lnN ')
        elif year == "2023postBPix":
            file.write('CMS_eff_m_trigger_2023BPix lnN ')
        else:
            file.write('CMS_eff_m_trigger_'+year+' lnN ')

        for i in range(nBins+4): # All except ZX
            file.write(trig_mu[year+'_'+channel]+' ')
        file.write('-\n')

    # Electrons efficiencies are decorrelated between RECO/ID, statistical and systematic components (syst are correlated between years)
    if channel == '4e' or channel == '2e2mu':
        
        if year == "2023preBPix":
            file.write('CMS_eff_e_id_2023 lnN ')
        elif year == "2023postBPix":
            file.write('CMS_eff_e_id_2023BPix lnN ')
        else:
            file.write('CMS_eff_e_id_'+year+' lnN ')
        for i in range(nBins+4):
            file.write(eff_e_id_stat[year+'_'+channel]+' ')
        file.write('-\n')
        
        if year == "2023preBPix":
            file.write('CMS_eff_e_reco_2023 lnN ')
        elif year == "2023postBPix":
            file.write('CMS_eff_e_reco_2023BPix lnN ')
        else:
            file.write('CMS_eff_e_reco_'+year+' lnN ')
        for i in range(nBins+4):
            file.write(eff_e_reco_stat[year+'_'+channel]+' ')
        file.write('-\n')

        file.write('CMS_eff_e_id lnN ')
        for i in range(nBins+4):
            file.write(eff_e_id_syst[year+'_'+channel]+' ')
        file.write('-\n')

        # --- RECO systematic (correlated across years) ---
        file.write('CMS_eff_e_reco_13p6TeV lnN ')
        for i in range(nBins+4):
            file.write(eff_e_reco_syst[year+'_'+channel]+' ')
        file.write('-\n')

        if year == "2023preBPix":
            file.write('CMS_eff_e_trigger_2023 lnN ')
        elif year == "2023postBPix":
            file.write('CMS_eff_e_trigger_2023BPix lnN ')
        else:
            file.write('CMS_eff_e_trigger_'+year+' lnN ')
        for i in range(nBins+4): # All except ZX
            file.write(trig_e[year+'_'+channel]+' ')
        file.write('-\n')

    # ZX
    file.write('CMS_HIG25015_hzz'+channel+'_Zjets_'+year+' lnN ')
    # for i in range(nBins+4): # All except ZX
    for i in range(nBins+4): # All except ZX
        file.write('- ')
    file.write(ZX[year+'_'+channel]+'\n')

    # Gaussian-constrained nuisance parameter , SMEARING systematic uncertainties
    if(channelNumber != 2):
        file.write('CMS_HIG25015_zz4l_mean_m_sig param 0.0 1.0\n')
        file.write('CMS_HIG25015_zz4l_sigma_m_sig param 0.0 0.03 [-1,1]\n') 
    if(channelNumber != 1):
        file.write('CMS_HIG25015_zz4l_mean_e_sig param 0.0 1.0\n')
        file.write('CMS_HIG25015_zz4l_sigma_e_sig param 0.0 0.1 [-1,1]\n')

    file.write('CMS_HIG25015_zz4l_n_sig_'+str(channelNumber)+'_'+year+' param 0.0 0.05\n')

    # Theoretical at 13.6 TeV taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap
    if not zzfloating:
        file.write('QCDscale_ggVV lnN ')
        for i in range(nBins+3): # Signal + out + fake + qqzz
            file.write('- ')
        file.write('1.039/0.961 -\n') #ggF (N3LO QCD + NLO EW), TH Gaussian % (+-3.9%)
        file.write('QCDscale_VV lnN ')
        for i in range(nBins+2): # Signal + out + fake
            file.write('- ')
        file.write('1.0325/0.958 - -\n')
        file.write('pdf_gg lnN ')
        for i in range(nBins+3): # Signal + out + fake + qqzz
            file.write('- ')
        file.write('1.032/0.968 -\n') #ggF (N3LO QCD + NLO EW), PDF+as% (+-3.2%)
        file.write('pdf_qqbar lnN ')
        for i in range(nBins+2): # Signal + out + fake
            file.write('- ')
        file.write('1.031/0.966 - -\n')
        file.write('CMS_HIG25015_kfactor_ggzz lnN ')
        for i in range(nBins+3): # Signal + out + fake  + bkg_qqzz
            file.write('- ')
        file.write('1.1 -\n')

    # JES
    if jes == True:

        if year == "2023preBPix":
            year = "2023"
        if year == "2023postBPix":
            year = "2023BPix"

        if obsName == "TCjmax": obsName_jes = "TCjMax"
        elif obsName == "TBjmax": obsName_jes = "TBjMax"
        elif obsName == "TCjmax_pT4l": obsName_jes = "TCjMax_ZZPt"
        else: obsName_jes = obsName

        if 'zzfloating' in obsName_jes: obsName_jes = obsName_jes.rsplit('_', 1)[0]

        for index,jesName in enumerate(jesNames_datacard):
            file.write('CMS_scale_j_'+jesName+' lnN ')
            for i in range(nBins+2): # Signals + out + fake
                file.write(str(fixJes(jesnp['signal_'+jesNames_datacard[index]+'_'+channel+'_'+year+'_'+obsName_jes.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)],
                                      jes_evts_noWeight['signal_'+jesNames_datacard[index]+'_'+channel+'_'+year+'_'+obsName_jes.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)])))
            file.write(str(fixJes(jesnp['qqzz_'+jesNames_datacard[index]+'_'+channel+'_'+year+'_'+obsName_jes.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)],
                                jes_evts_noWeight['qqzz_'+jesNames_datacard[index]+'_'+channel+'_'+year+'_'+obsName_jes.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)])))
            file.write(str(fixJes(jesnp['ggzz_'+jesNames_datacard[index]+'_'+channel+'_'+year+'_'+obsName_jes.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)],
                                jes_evts_noWeight['ggzz_'+jesNames_datacard[index]+'_'+channel+'_'+year+'_'+obsName_jes.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)])))
            #file.write(str(fixJes(jesnp['ZX_'+jesNames[index]+'_'+channel+'_'+year+'_'+obsName.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)],
            #                      jes_evts_noWeight['ZX_'+jesNames[index]+'_'+channel+'_'+year+'_'+obsName.replace('pT4l', 'ZZPt')+'_recobin'+str(obsBin)]))+'\n')

            file.write('- \n') # none for Z+X


            # file.write('CMS_scale_j_ZX lnN ')
            # for i in range(nBins+4): # All except ZX
            #     file.write('- ')
            # file.write(str(jesnp['ZX_'+channel+'_'+obsName+'_genbin'+str(obsBin)+'_recobin'+str(obsBin)])+'\n')

                #file.write('JES param 0.0 1.0\n')

    file.close()

    print(os.getcwd())

#This createDatacard_ggH function is obsolete!! It's not used anymore
def createDatacard_ggH(obsName, channel, nBins, obsBin, observableBins, physicalModel, year, nData, jes, lowerBound, upperBound, yearSetting):
    # Name of the bin (aFINALSTATE_ recobinX)
    if(channel == '4mu'): channelNumber = 1
    if(channel == '4e'): channelNumber = 2
    if(channel == '2e2mu'): channelNumber = 3
    # binName = 'a'+str(channelNumber)+'_recobin'+str(obsBin)

    # _recobin = str(observableBins[obsBin]).replace('.', 'p') + '_' + str(observableBins[obsBin+1]).replace('.', 'p')
    _recobin = str(observableBins[obsBin]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[obsBin+1]).replace('.', 'p').replace('-','m')
    if int(observableBins[obsBin+1]) > 1000:
        _recobin = 'GT'+str(int(observableBins[obsBin]))

    _obsName = {'pT4l': 'PTH', 'rapidity4l': 'YH', 'pTj1': 'pTj1', 'Nj': 'Nj'}
    binName = 'hzz_' + _obsName[obsName] + '_' + _recobin + '_cat' + channel
    # Root of the name of the process (signal from genBin)
    # processName = 'smH'+channel+'Bin'
    processName = 'ggH_' + _obsName[obsName]
    xHName = 'xH_' + _obsName[obsName]
    # Background expectations in [105,160]
    sys.path.append('../inputs')
    #_temp = __import__('inputs_bkgTemplate_'+obsName, globals(), locals(), ['expected_yield'], -1)
    _temp = __import__('inputs_bkgTemplate_'+obsName, globals(), locals(), ['expected_yield'], 0) # spencer
    expected_yield = _temp.expected_yield
    if jes:

        sys.path.append('../coefficients/JES')
        #_temp = __import__('JESNP_'+obsName+'_'+str(year), globals(), locals(), ['JESNP'], -1)
        _temp = __import__('JESNP_'+obsName+'_'+str(year), globals(), locals(), ['JESNP'], 0) # spencer
        jesnp = _temp.JESNP

        if year == "2023preBPix":
            year_jes = "2023"
        elif year == "2023postBPix":
            year_jes = "2023BPix"
        else:
            year_jes = year 

        jesNames = ['Absolute', 'Absolute_year', 'BBEC1', 'BBEC1_year', 'EC2', 'EC2_year', 'FlavorQCD', 'HF', 'HF_year', 'RelativeBal', 'RelativeSample_year']
        jesNames_datacard = [j.replace('year',year_jes) for j in jesNames] # The name of the nuisance in the datacard should have the correspoding year


        sys.path.remove('../coefficients/JES')
        # print(jesnp)


    sys.path.remove('../inputs')

    # lumi
    if yearSetting == 'Full':
        lumi = {}
        lumi['2016'] = '1.01'
        lumi['2017'] = '1.02'
        lumi['2018'] = '1.015'
        lumi_corr_16_17_18 = {}
        lumi_corr_16_17_18['2016'] = '1.006'
        lumi_corr_16_17_18['2017'] = '1.009'
        lumi_corr_16_17_18['2018'] = '1.02'
        lumi_corr_17_18 = {}
        lumi_corr_17_18['2017'] = '1.006'
        lumi_corr_17_18['2018'] = '1.002'
    elif yearSetting == 'Run3':
        # RUN III correlated luminosity scheme (2022/23/24)
        lumi = {
        '2022':'1.0138',
        '2022EE':'1.0138',
        '2023preBPix':'1.0017',
        '2023postBPix':'1.0017',
        '2024':'1.0020',
        }
        # Additional reduced nuisances for Run-3:
        # lumi_13p6TeV_2324 applies to 2023 & 2024 / lumi_13p6TeV_2024 applies to 2024 only
        lumi_extra = {
        'lumi_13p6TeV_2324': {'2023preBPix': '1.0127', '2023postBPix': '1.0127', '2024': '1.0068'},
        'lumi_13p6TeV_2024': {'2024': '1.0144'}
        }
    else:
        lumi = {}
        lumi['2016'] = '1.026'
        lumi['2017'] = '1.025'
        lumi['2018'] = '1.023'
        lumi['2022'] = '1.014'
        lumi['2022EE'] = '1.014'
        lumi['2023preBPix'] = '1.013'
        lumi['2023postBPix'] = '1.013'
        lumi['2024'] = '1.013' # using 2023postBPix value as placeholder
        
    # Lepton efficiency
    # Values taken from:
    # https://indico.cern.ch/event/1125278/contributions/4723319/attachments/2420259/4142589/LepSyst.pdf
    # Computed with:
    # https://github.com/asculac/LeptonSystematics
    eff_mu = {}
    eff_mu['2016_2e2mu'] = '0.986/1.007'
    eff_mu['2016_4mu'] = '0.981/1.01'
    eff_mu['2017_2e2mu'] = '0.986/1.006'
    eff_mu['2017_4mu'] = '0.981/1.009'
    eff_mu['2018_2e2mu'] = '0.986/1.006'
    eff_mu['2018_4mu'] = '0.981/1.008'

    #TODO: Update 2023 eff_mu values. The values for 2023 preBPix/postBPix are copy-pasted from 2022 at the moment
    eff_mu['2023preBPix_2e2mu'] = '0.986/1.008' #spencer
    eff_mu['2023preBPix_4mu'] = '0.981/1.01' #spencer
    eff_mu['2023postBPix_2e2mu'] = '0.986/1.008' # spencer
    eff_mu['2023postBPix_4mu'] = '0.981/1.01' # spencer

    eff_mu['2023preBPix_2e2mu'] = '0.987/1.005' # ok
    eff_mu['2023preBPix_4mu'] = '0.982/1.006' # ok
    eff_mu['2023postBPix_2e2mu'] = '0.987/1.006' # ok
    eff_mu['2023postBPix_4mu'] = '0.982/1.007' # ok

    eff_mu['2024_2e2mu'] = '0.988/1.003' # ok
    eff_mu['2024_4mu'] = '0.996/1.017' # ok
    
    eff_e = {}
    eff_e['2016_2e2mu'] = '0.934/1.062'
    eff_e['2016_4e'] = '0.891/1.093'
    eff_e['2017_2e2mu'] = '0.953/1.043'
    eff_e['2017_4e'] = '0.915/1.064'
    eff_e['2018_2e2mu'] = '0.95/1.052'
    eff_e['2018_4e'] = '0.905/1.077'

    #TODO: Update 2023 eff_e values. The values for 2023 preBPix/postBPix are copy-pasted from 2022 at the moment
    eff_e['2023preBPix_2e2mu'] = '0.927/1.069' # spencer
    eff_e['2023preBPix_4e'] = '0.884/1.103' # spencer
    eff_e['2023postBPix_2e2mu'] = '0.927/1.069' # spencer
    eff_e['2023postBPix_4e'] = '0.884/1.103' # spencer

    eff_e['2023preBPix_2e2mu'] = '0.882/1.116' # ok
    eff_e['2023preBPix_4e'] = '0.814/1.177' # ok
    eff_e['2023postBPix_2e2mu'] = '0.852/1.146' # ok
    eff_e['2023postBPix_4e'] = '0.778/1.215' # ok

    eff_e['2024_2e2mu'] = '0.933/1.064'# ok
    eff_e['2024_4e'] = '0.887/1.098'# ok

    
    # ZX
    # Values taken from:
    # https://indico.cern.ch/event/1109103/contributions/4799990/attachments/2415431/4133107/HIG21009_HZZreport.pdf
    # Computed with:
    # https://github.com/CJLST/ZZAnalysis/blob/Run2_CutBased_UL/AnalysisStep/test/ZpXEstimation/
    ZX = {}
    ZX['2016_2e2mu'] = '0.756295/1.25114'
    ZX['2016_4e'] = '0.598182/1.43059'
    ZX['2016_4mu'] = '0.694678/1.30555'
    ZX['2017_2e2mu'] = '0.765736/1.236385'
    ZX['2017_4e'] = '0.646521/1.36398'
    ZX['2017_4mu'] = '0.694063/1.30623'
    ZX['2018_2e2mu'] = '0.7660052/1.235647'
    ZX['2018_4e'] = '0.650486/1.35893'
    ZX['2018_4mu'] = '0.69554/1.30465'

    #TODO: 2024 is evaluated using 2023 DY, should be recomputed
    ZX['2023preBPix_2e2mu'] = '0.724/1.263' # spencer
    ZX['2023preBPix_4e'] = '0.575/1.398' # spencer
    ZX['2023preBPix_4mu'] = '0.677/1.321' # spencer
    ZX['2023postBPix_2e2mu'] = '0.724/1.263' # spencer
    ZX['2023postBPix_4e'] = '0.575/1.398' # spencer
    ZX['2023postBPix_4mu'] = '0.677/1.321' # spencer

    ZX['2024_2e2mu'] = '0.766/1.232' # spencer 
    ZX['2024_4e'] = '0.654/1.339' # spencer 
    ZX['2024_4mu'] = '0.696/1.303' # spencer
    
    # -------------------------------------------------------------------------------------------------

    file = open('../datacard/datacard_'+year+'/hzz4l_GGH_'+channel+'S_13TeV_xs_'+_obsName[obsName]+'_bin'+str(obsBin)+'_'+physicalModel+'.txt', 'w+')

    file.write('imax 1 \n')
    file.write('jmax * \n')
    file.write('kmax * \n')

    file.write('------------ \n')

    file.write('shapes * * hzz4l_'+channel+'S_13TeV_xs_SM_125_'+_obsName[obsName]+'_'+physicalModel+'.Databin'+str(obsBin)+'.root w:$PROCESS\n')

    file.write('------------ \n')

    file.write('bin '+binName+'\n')
    file.write('observation '+str(nData)+'\n')

    file.write('------------ \n')
    file.write('## mass window ['+str(lowerBound)+','+str(upperBound)+']\n')
    file.write('bin ')
    # for i in range(nBins+5): # In addition to the observableBins, there are OutsideAcceptance, fakeH, bkg_ggzz, bkg_qqzz, bkg_zjets
    for i in range(2*nBins+5):
        file.write(binName+' ')
    file.write('\n')
    file.write('process ')
    if '_' in obsName and not 'floating' in obsName and not 'kL' in obsName and not obsName == 'Nj':
        for i in range(nBins):
            file.write(processName+'_'+str(observableBins[i][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][1]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][2]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][3]).replace('.', 'p').replace('-','m')+' ')
        for i in range(nBins):
            file.write(xHName+'_'+str(observableBins[i][0]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i][1]).replace('.', 'p').replace('-','m')+str(observableBins[i][2]).replace('.', 'p').replace('-','m')+str(observableBins[i][3]).replace('.', 'p').replace('-','m')+' ')
    else:
        for i in range(nBins):
            # file.write(processName+str(i)+' ')
            if observableBins[i+1] > 1000:
                file.write(processName+'_GT'+str(int(observableBins[i]))+' ')
            else:
                file.write(processName+'_'+str(observableBins[i]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i+1]).replace('.', 'p').replace('-','m')+' ')
        for i in range(nBins):
            # file.write(processName+str(i)+' ')
            if observableBins[i+1] > 1000:
                file.write(xHName+'_GT'+str(int(observableBins[i]))+' ')
            else:
                file.write(xHName+'_'+str(observableBins[i]).replace('.', 'p').replace('-','m')+'_'+str(observableBins[i+1]).replace('.', 'p').replace('-','m')+' ')
    file.write('OutsideAcceptance nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    #file.write('nonResH bkg_qqzz bkg_ggzz bkg_zjets')
    file.write('\n')
    file.write('process ')
    for i in range(nBins):
        file.write('-'+str(i+1)+' ')
    # nonResH bkg_qqzz bkg_ggzz bkg_zjets
    file.write('1 2 3 4 5')
    # xH as bkg processBin
    for i in range(nBins):
        file.write(str(i+5)+' ')
    file.write('\n')
    file.write('rate ')
    # ggH, XH bins and nonRes contribution
    for i in range(2*nBins+2):
        file.write('1.0 ')
    # file.write(str(bkg_qqzz[year+'_'+channel])+' '+str(bkg_ggzz[year+'_'+channel])+' '+str(bkg_zx[year+'_'+channel])+'\n') #Old implementation with hard coding bkg expectation values
    file.write(str(expected_yield[int(year),'qqzz',channel])+' '+str(expected_yield[int(year),'ggzz',channel])+' '+str(expected_yield[int(year),'ZX',channel])+'\n')
    file.write('------------ \n')

    # norm_fake
    file.write('norm_nonResH lnU ')
    # ggH + xH
    for i in range(2*nBins+1):
        file.write('- ')
    file.write('10.0 - - -    # [/10,*10]\n')

    if yearSetting == 'Full':
        # lumi_uncorrelated
        file.write('lumi_13TeV_'+year+'_uncorrelated lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+4): # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n') # ZX
        # lumi_correlated_16_17_18
        file.write('lumi_13TeV_correlated_16_17_18 lnN ')
        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+4): # All except ZX
            file.write(lumi_corr_16_17_18[year]+' ')
        file.write('-\n') # ZX
        # lumi_correlated_17_18
        if year == '2017' or year == '2018':
            file.write('lumi_13TeV_correlated_17_18 lnN ')
            # for i in range(nBins+4): # All except ZX
            for i in range(2*nBins+4): # All except ZX
                file.write(lumi_corr_17_18[year]+' ')
            file.write('-\n') # ZX
            # lumi_correlated_16_17_18
            file.write('lumi_13TeV_correlated_16_17_18 lnN ')
            for i in range(nBins+4): # All except ZX
                file.write(lumi_corr_16_17_18[year]+' ')
            file.write('-\n') # ZX
            # lumi_correlated_17_18
            if year == '2017' or year == '2018':
                file.write('lumi_13TeV_correlated_17_18 lnN ')
                for i in range(nBins+4): # All except ZX
                    file.write(lumi_corr_17_18[year]+' ')
                file.write('-\n') # ZX

    elif yearSetting == 'Run3':

        # correlated Run-3 lumi nuisances
        # lumi_13p6TeV_222324 (2022,2023,2024)
        file.write('lumi_13p6TeV_222324 lnN ')
        for i in range(2*nBins+4): 
            # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n')

        # lumi_13p6TeV_2324 if applicable
        if year in lumi_extra['lumi_13p6TeV_2324']:
            file.write('lumi_13p6TeV_2324 lnN ')
            for i in range(2*nBins+4):
                file.write(lumi_extra['lumi_13p6TeV_2324'][year]+' ')
            file.write('-\n')

        # lumi_13p6TeV_2024 if applicable
        if year in lumi_extra['lumi_13p6TeV_2024']:
            file.write('lumi_13p6TeV_2024 lnN ')
            for i in range(2*nBins+4):
                file.write(lumi_extra['lumi_13p6TeV_2024'][year]+' ')
            file.write('-\n')

    else:

        # lumi
        if year == '2022' or year == '2022EE':
            file.write('lumi_13p6TeV_2022 lnN ')
        elif year == '2023preBPix' or year == '2023postBPix':
            file.write('lumi_13p6TeV_2023 lnN ')
        elif year == '2024':
            file.write('lumi_13p6TeV_2024 lnN ')
        else:
            file.write('lumi_13TeV_'+year+' lnN ')

        # for i in range(nBins+4): # All except ZX
        for i in range(2*nBins+4): # All except ZX
            file.write(lumi[year]+' ')
        file.write('-\n') # ZX

    # Lepton efficiency
    if channel == '4mu' or channel == '2e2mu':
        file.write('CMS_eff_m lnN ')
        # for i in range(nBins+4)