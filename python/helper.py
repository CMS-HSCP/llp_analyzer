import numpy as np
import json
import sys
sys.path.append('/storage/af/user/christiw/gpu/christiw/llp/delayed_jet_analyzer/lib/')
from histo_utilities import create_TGraph



path = '/storage/af/user/christiw/login-1/christiw/LLP/dedx/CMSSW_10_6_30/src/llp_analyzer/data/xsec/json/'
filenames = {
   'gluino': 'pp13_gluino_NNLO+NNLL.json',
#     'pp13_slep_LR_NLO+NLL_PDF4LHC.json',
#     'stop':'pp13_squark_NNLO+NNLL.json',
    'chargino':'pp13_wino_C1C1_NLO+NLL.json',
    'stau': 'pp13_stau_LR_NLO+NLL_PDF4LHC.json',
}
mass = {}
theoretical_xsec = {}
for f, file in filenames.items():
    data = json.load(open( path + file))
    mass[f] = []
    theoretical_xsec[f] = []
    for k,v in data['data'].items():
        mass[f].append(int(k))
        theoretical_xsec[f].append(float(v['xsec_pb']))
    mass[f] = np.array(mass[f])
    theoretical_xsec[f] = np.array(theoretical_xsec[f])
    inds = mass[f].argsort()
    mass[f] = mass[f][inds]
    theoretical_xsec[f] = theoretical_xsec[f][inds]
def get_theoretical_xsec(m, signal):
    if signal == 'stau':
        h = create_TGraph(mass[signal],theoretical_xsec[signal])
        return h.Eval(m)
    elif signal == 'gluino':
        h = create_TGraph(mass[signal],theoretical_xsec[signal])
        return h.Eval(m)
    else:
        assert(False)

def find_intersect(graph1,graph2):
    
    mass = list(np.arange(0,2500))
    diff = []
    for m in mass:
        diff.append(graph1.Eval(m)-graph2.Eval(m))
    diff = np.sign(np.array(diff))
    for i in range(len(diff)):
        if not diff[i] == diff[i+1]:
            return mass[i]




def make_datacard_2tag(outDataCardsDir,modelName,  signal_rate, bkg_rate, observation, sig_unc, norm, prefix):
    a,b,c = bkg_rate[0], bkg_rate[1], bkg_rate[2]
    nSig = len(signal_rate.keys())
    nBins = 3
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# norm {0} \n'.format(norm))


    text_file.write('imax {0} \n'.format(nBins))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')



    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f}  \n'.format(observation[0],observation[1],observation[2]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * nBins + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * nBins + '\n')
    rate_string = 'rate'
    for i in range(nBins):# 3 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')

    text_file.write(prefix+'A    rateParam       chA     bkg      (@0/(2*@1))**2*@1                     '+prefix+'B,'+prefix+'C \n')
    text_file.write(prefix+'B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7))
    text_file.write(prefix+'C   rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7))

    for k,v in signal_rate.items():
        text_file.write('norm rateParam * {0} 1  \n'.format(k))
#   #### uncertainties ####
    for k,v in sig_unc.items():
        unc_text = k+' \t lnN'
        for ele in v:
            if ele == 0.0:unc_text += ' \t -'
            else: unc_text += ' \t '+str(ele+1)
            unc_text += '\t - '
        
        text_file.write(unc_text + ' \n')
#     for i in range(len(bkg_unc_name)):
#         bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(4*nSig+3) + '\t ' + str(1+bkg_unc[i]) + ' \n'
#         text_file.write(bkg_unc_text)

    text_file.close()

def make_datacard(outDataCardsDir,modelName,  signal_rate, bkg_rate, observation, sig_unc, norm, prefix):
    a,b,c,d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    c1 = a/b
    c2 = c/b
    nSig = len(signal_rate.keys())
    text_file = open(outDataCardsDir+modelName+".txt", "w")
    text_file.write('# norm {0} \n'.format(norm))


    text_file.write('imax {0} \n'.format(4))
    text_file.write('jmax {0} \n'.format(nSig))
    text_file.write('kmax * \n')
    text_file.write('shapes * * FAKE \n')



    text_file.write('--------------- \n')
    text_file.write('--------------- \n')
    text_file.write('bin \t chA \t chB \t chC \t chD \n')
    text_file.write('observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n'.format(observation[0],observation[1],observation[2],observation[3]))
    text_file.write('------------------------------ \n')
    text_file.write('bin '+'\t chA ' * (1+nSig) + '\t chB ' * (1+nSig) +'\t chC '*(1+nSig) +'\t chD '*(1+nSig) +'\n')
    process_name = '\t '+ (' \t ').join(list(signal_rate.keys())) + '\t bkg '
    text_file.write('process ' + process_name * 4 + '\n')
    process_number = '\t '+ (' \t ').join(list((np.arange(nSig)*-1).astype(str))) + '\t 1'
    text_file.write('process ' + process_number * 4 + '\n')
    rate_string = 'rate'
    for i in range(4):# 4 bins
        for k,v in signal_rate.items():
            rate_string +='\t {0:e} '.format(v[i])
        rate_string += '\t 1 '
    text_file.write(rate_string+'\n')
    text_file.write('------------------------------ \n')
    text_file.write(prefix+'A    rateParam       chA     bkg      (@0*@2/@1)                    '+prefix+'B,'+prefix+'C,'+prefix+'D \n')
    text_file.write(prefix+'B    rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(b, b*7))
    text_file.write(prefix+'C    rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(c, c*7))
    text_file.write(prefix+'D    rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n'.format(d, d*7))
    
    for k,v in signal_rate.items():
        text_file.write('norm rateParam * {0} 1  \n'.format(k))
#   #### uncertainties ####
    for k,v in sig_unc.items():
        unc_text = k+' \t lnN'
        for ele in v:
            if ele == 0.0:unc_text += ' \t -'
            else: unc_text += ' \t '+str(ele+1)
            unc_text += '\t - '
        
        text_file.write(unc_text + ' \n')
#     for i in range(len(bkg_unc_name)):
#         bkg_unc_text = bkg_unc_name[i] + ' \t lnN ' + '\t - '*(4*nSig+3) + '\t ' + str(1+bkg_unc[i]) + ' \n'
#         text_file.write(bkg_unc_text)

    text_file.close()
