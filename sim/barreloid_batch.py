"""
barreloid_batch.py

MPI usage:
    mpiexec -n 8 nrniv -python -mpi barreloid_batch.py
    mpiexec -n 64 nrniv -python -mpi barreloid_batch.py

Contributors: joao.moreira@downstate.edu, salvadordura@gmail.com
"""
import  os
import  sys 
import  numpy           as      np
from    netpyne.batch   import  Batch
import  NetPyNE_BBP
import  json

def batchWeight(batch_label = None, session_type = 'detailed', batch_flag=''):
    if batch_label is None: 
        print(' --- Select one of the possible simulation types:')
        print(  
                '\n\t_________________________________________________',
                '\n\tbarreloid_batch_fig02_sleep_wake',
                '\n\t_________________________________________________',
                '\n\tbarreloid_batch_fig05_awake',
                '\n\t_________________________________________________',
                '\n\tbarreloid_batch_fig06_sleep',
                '\n\t_________________________________________________',
                '\n\tbarreloid_batch_fig07e_0_uniform_ct_modulation_0\n\tbarreloid_batch_fig07e_0_uniform_ct_modulation_10\n\tbarreloid_batch_fig07e_0_uniform_ct_modulation_25\n\tbarreloid_batch_fig07e_0_uniform_ct_modulation_50\n\tbarreloid_batch_fig07e_0_uniform_ct_modulation_100',
                '\n\tbarreloid_batch_fig07e_1_mixed_ct_modulation_0\n\tbarreloid_batch_fig07e_1_mixed_ct_modulation_10\n\tbarreloid_batch_fig07e_1_mixed_ct_modulation_25\n\tbarreloid_batch_fig07e_1_mixed_ct_modulation_50\n\tbarreloid_batch_fig07e_1_mixed_ct_modulation_100',
                '\n\tbarreloid_batch_fig07e_2_closed_ct_modulation_0\n\tbarreloid_batch_fig07e_2_closed_ct_modulation_10\n\tbarreloid_batch_fig07e_2_closed_ct_modulation_25\n\tbarreloid_batch_fig07e_2_closed_ct_modulation_50\n\tbarreloid_batch_fig07e_2_closed_ct_modulation_100',
                '\n\t_________________________________________________\n\n',
                '\n\tbarreloid_batch_fig07f_0_uniform_ct_modulation_0\n\tbarreloid_batch_fig07f_0_uniform_ct_modulation_10\n\tbarreloid_batch_fig07f_0_uniform_ct_modulation_25\n\tbarreloid_batch_fig07f_0_uniform_ct_modulation_50\n\tbarreloid_batch_fig07f_0_uniform_ct_modulation_100',
                '\n\tbarreloid_batch_fig07f_1_mixed_ct_modulation_0\n\tbarreloid_batch_fig07f_1_mixed_ct_modulation_10\n\tbarreloid_batch_fig07f_1_mixed_ct_modulation_25\n\tbarreloid_batch_fig07f_1_mixed_ct_modulation_50\n\tbarreloid_batch_fig07f_1_mixed_ct_modulation_100',
                '\n\tbarreloid_batch_fig07f_2_closed_ct_modulation_0\n\tbarreloid_batch_fig07f_2_closed_ct_modulation_10\n\tbarreloid_batch_fig07f_2_closed_ct_modulation_25\n\tbarreloid_batch_fig07f_2_closed_ct_modulation_50\n\tbarreloid_batch_fig07f_2_closed_ct_modulation_100',
                '\n\t_________________________________________________\n\n',
                )
        return

    ##################################################################################################################################################################

    #############################################
    # Brainstem Inputs - Hartings               #
    #############################################

    # 25 uS - even deflection times (all = 250 ms)
    bs_wake_none    = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|None_samplingBins|25us_stims|None_simplifiedDataset.pkl'
    bs_wake_0       = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|0_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_45      = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|45_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_90      = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|90_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_135     = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|135_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_180     = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|180_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_225     = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|225_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_270     = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|270_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'
    bs_wake_315     = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/mleSpikes71_deflectionAngle|315_samplingBins|25us_stims|2000@2250|3000@3250|4000@4250|5000@5250|6000@6250|7000@7250|8000@8250|9000@9250|10000@10020_simplifiedDataset.pkl'

    bs_sleep_none   = '../stims/Artificial/CT/spiking_dict__nCells_1000__sampling_25_uS__freq_10_Hz__downTimes__None__.json'

    ##########################
    # Corticothalamic Inputs #
    ##########################
    ct_constant={
                    'uniform':      '../stims/Artificial/CT/v3_spiking_dict__nCells_818__sampling_25_uS__freq_5_Hz__downTimes__None__.json',
                }
    # CT up-down
    ct_upDown2={
                    'uniform':      '../stims/Artificial/CT/v4_spiking_dict__nCells_818__sampling_25_uS__freq_5_Hz__downFreq_1__downTimes__0|2500_3000|4500_5000|6500_7000|8500_9000|10000__.json',
                }
    
    ###########################
    # Hartings Angular Tuning #
    ###########################
    # 10%   -   Replace with 82 angular tuned CT cells - Hartings Angular Tuning
    replace_ct__10pct__None = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|None_samplingBins|25us_stims|None_simplifiedDataset.pkl'
    replace_ct__10pct__0    = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|0_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__45   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|10_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__90   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|20_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__135  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|31_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__180  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|41_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__225  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|51_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__270  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|62_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__10pct__315  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes60_deflectionAngle|72_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'

    # 25%   -   Replace with 204 angular tuned CT cells - Hartings Angular Tuning
    replace_ct__25pct__None = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|None_samplingBins|25us_stims|None_simplifiedDataset.pkl'
    replace_ct__25pct__0    = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|0_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__45   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|26_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__90   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|51_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__135  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|76_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__180  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|102_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__225  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|128_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__270  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|153_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__25pct__315  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes61_deflectionAngle|178_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'

    # 50%   -   Replace with 409 angular tuned CT cells - Hartings Angular Tuning
    replace_ct__50pct__None = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|None_samplingBins|25us_stims|None_simplifiedDataset.pkl'
    replace_ct__50pct__0    = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|0_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__45   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|51_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__90   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|102_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__135  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|153_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__180  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|204_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__225  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|256_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__270  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|307_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__50pct__315  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes62_deflectionAngle|358_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'

    # 100%   -   Replace with 818 angular tuned CT cells - Hartings Angular Tuning
    replace_ct__100pct__None = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|None_samplingBins|25us_stims|None_simplifiedDataset.pkl'
    replace_ct__100pct__0    = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|0_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__45   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|102_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__90   = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|204_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__135  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|307_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__180  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|409_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__225  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|511_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__270  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|614_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'
    replace_ct__100pct__315  = '../stims/Minnery_2003/fig1/deflection_model/spike_dicts/ctSpikes63_deflectionAngle|716_samplingBins|25us_stims|2008.73@2258.73|3008.73@3258.73|4008.73@4258.73|5008.73@5258.73|6008.73@6258.73|7008.73@7258.73|8008.73@8258.73|9008.73@9258.73|10000@10020_simplifiedDataset.pkl'

    ############################################################################################################################################################################
    # === Paper simulations === #
    
    # === Sleep-Wake === #
    # --- barreloid_batch_fig02_sleep_wake
    if      batch_label == 'barreloid_batch_fig02_sleep_wake':
        network_state_dict={
            'sleep__baseline_5':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 3.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__baseline_10Hz':    {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #        
        }
    
    # === Awake === #
    # --- barreloid_batch_fig05_awake
    elif    batch_label == 'barreloid_batch_fig05_awake':
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
        }
    
    # === Sleep === #
    # --- barreloid_batch_fig06_sleep
    elif    batch_label == 'barreloid_batch_fig06_sleep':
        network_state_dict={
                'sleep__baseline_5':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 3.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_upDown2, 'replace_ct_dict': None    },  #
        }

    # === Awake CT modulation === #
    # --- barreloid_batch_fig07e_0_uniform_ct_modulation_0 / barreloid_batch_fig07e_1_mixed_ct_modulation_0 / barreloid_batch_fig07e_2_closed_ct_modulation_0
    elif    (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_0'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
        }

    # --- barreloid_batch_fig07e_0_uniform_ct_modulation_10 / barreloid_batch_fig07e_1_mixed_ct_modulation_10 / barreloid_batch_fig07e_2_closed_ct_modulation_10
    elif    (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_10'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__315    },  #
        }

    # --- barreloid_batch_fig07e_0_uniform_ct_modulation_25 / barreloid_batch_fig07e_1_mixed_ct_modulation_25 / barreloid_batch_fig07e_2_closed_ct_modulation_25
    elif    (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_25'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__315    },  #
        }

    # --- barreloid_batch_fig07e_0_uniform_ct_modulation_50 / barreloid_batch_fig07e_1_mixed_ct_modulation_50 / barreloid_batch_fig07e_2_closed_ct_modulation_50
    elif    (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_50'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__315    },  #
        }
    
    # --- barreloid_batch_fig07e_0_uniform_ct_modulation_100 / barreloid_batch_fig07e_1_mixed_ct_modulation_100 / barreloid_batch_fig07e_2_closed_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_100') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_100') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_100'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__315    },  #
        }

    # === Awake CT modulation === #
    # --- barreloid_batch_fig07f_0_uniform_ct_modulation_0 / barreloid_batch_fig07f_1_mixed_ct_modulation_0 / barreloid_batch_fig07f_2_closed_ct_modulation_0
    elif    (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_0'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
        }

    # --- barreloid_batch_fig07f_0_uniform_ct_modulation_10 / barreloid_batch_fig07f_1_mixed_ct_modulation_10 / barreloid_batch_fig07f_2_closed_ct_modulation_10
    elif    (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_10'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__10pct__315    },  #
        }

    # --- barreloid_batch_fig07f_0_uniform_ct_modulation_25 / barreloid_batch_fig07f_1_mixed_ct_modulation_25 / barreloid_batch_fig07f_2_closed_ct_modulation_25
    elif    (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_25'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__25pct__315    },  #
        }

    # --- barreloid_batch_fig07f_0_uniform_ct_modulation_50 / barreloid_batch_fig07f_1_mixed_ct_modulation_50 / barreloid_batch_fig07f_2_closed_ct_modulation_50
    elif    (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_50'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__50pct__315    },  #
        }
    
    # --- barreloid_batch_fig07f_0_uniform_ct_modulation_100 / barreloid_batch_fig07f_1_mixed_ct_modulation_100 / barreloid_batch_fig07f_2_closed_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_100') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_100') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_100'):
        network_state_dict={
            'wake__0deg':             {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_0,     'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__0    },  #
            'wake__45deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_45,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__45    },  #
            'wake__90deg':            {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_90,    'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__90    },  #
            'wake__135deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_135,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__135    },  #
            'wake__180deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_180,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__180    },  #
            'wake__225deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_225,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__225    },  #
            'wake__270deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_270,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__270    },  #
            'wake__315deg':           {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_wake_315,   'ct_dict': ct_constant, 'replace_ct_dict': replace_ct__100pct__315    },  #
        }

    # # --- test_TRN_RMP
    # elif    batch_label == 'test_TRN_RMP':
    #     # network_state_dict={
    #     #     'sleep__baseline':        {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 6.0,  6.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            
    #     #     'sleep__baseline_1':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 4.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
    #     #     'sleep__baseline_2':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 4.0, 10.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
    #     #     'sleep__baseline_3':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 2.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
    #     #     'sleep__baseline_4':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 2.0, 10.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
    #     #     'sleep__baseline_5':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
    #     #     'sleep__baseline_6':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, 10.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #
            
    #     #     'wake__baseline_10Hz':    {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 0.5,   'rescale_GABA_g_extra': 0.3,  'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0, -5.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_constant, 'replace_ct_dict': None    },  #        
    #     # }
    #     # network_state_dict={
    #     #     'sleep__baseline_5':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 0.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_upDown2, 'replace_ct_dict': None    },  #
    #     # }
    #     # # 3_1
    #     # network_state_dict={
    #     #     'sleep__baseline_5':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 3.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_upDown2, 'replace_ct_dict': None    },  #
    #     # }
    #     # 3_2
    #     network_state_dict={
    #         'sleep__baseline_5':      {'w_CT_VPM':1.0 ,'w_CT_TRN':1.5, 'rescale_GABA_g_intra': 1.00,  'rescale_GABA_g_extra': 1.00, 'NMDA_L6_VPM': 2.0, 'NMDA_L6_TRN': 1.0, 'w_MLe_VPM': 3.0 , 'w_VPM_TRN':0.5, 'shift_RPM':[ 6.0,  8.0], 'bs_dict': bs_sleep_none, 'ct_dict': ct_upDown2, 'replace_ct_dict': None    },  #
    #     }

    ###########################################################################################################################################################################

    for state in network_state_dict.keys():
        if 'rescale_KLeak_TC_g' not in network_state_dict[state].keys(): network_state_dict[state].update({'rescale_KLeak_TC_g':0.0})
        if 'rescale_KLeak_RE_g' not in network_state_dict[state].keys(): network_state_dict[state].update({'rescale_KLeak_RE_g':0.0})

    for state in network_state_dict.keys():
        if 'ct_downTimes' not in network_state_dict[state].keys(): network_state_dict[state].update({'ct_downTimes':None})

    batch_groups = list(network_state_dict.keys())
    for key in batch_groups:
        if '__noGap' in key:    network_state_dict[key].update({'gap_gmax':0})
        else:                   network_state_dict[key].update({'gap_gmax':0.2e-6})

        # sleep-wake cell parameters
        if   'sleep__' in key:   network_state_dict[key].update({'rescale_iNap':1.00, 'rescale_ih':1.00, 'rescale_iT':1.00, 'rescale_iA':1.00})
        elif 'wake__'  in key:   network_state_dict[key].update({'rescale_iNap':0.25, 'rescale_ih':0.25, 'rescale_iT':0.25, 'rescale_iA':1.50})
        else:                    network_state_dict[key].update({'rescale_iNap':1.00, 'rescale_ih':1.00, 'rescale_iT':1.00, 'rescale_iA':1.00})


    ############################################################################################################################################################################
    parameter_space_rescale_inap_gmax = [network_state_dict[state]['rescale_iNap'] for state in network_state_dict.keys()]
    parameter_space_rescale_ih_gmax   = [network_state_dict[state]['rescale_ih']   for state in network_state_dict.keys()]
    parameter_space_rescale_it_gmax   = [network_state_dict[state]['rescale_iT']   for state in network_state_dict.keys()]
    parameter_space_rescale_ia_gmax   = [network_state_dict[state]['rescale_iA']   for state in network_state_dict.keys()]

    ############################################################################################################################################################################
    
    parameter_space_replace_spk_file = [network_state_dict[state]['replace_ct_dict']   for state in network_state_dict.keys()]
    
    ############################################################################################################################################################################

    parameter_space_rescale_iH_gmax_VPM = [1.0] # (default: 1.0)
    parameter_space_rescale_iH_gmax_TRN = [1.0] # (default: 1.0)

    ############################################################################################################################################################################

    parameter_space_addUSE_MLe      = [0]
    parameter_space_rescaleDep_MLe  = [0.75]

    ##################################################################################################################################################################

    parameter_space_KLeak_VPM    = [network_state_dict[state]['rescale_KLeak_TC_g']   for state in network_state_dict.keys()]
    parameter_space_KLeak_TRN    = [network_state_dict[state]['rescale_KLeak_RE_g']   for state in network_state_dict.keys()]
    parameter_space_TRN_VPM      = [network_state_dict[state]['rescale_GABA_g_extra'] for state in network_state_dict.keys()] # (Baghdoyan, 2012) reference for GABA rescaling
    parameter_space_TRN_TRN      = [network_state_dict[state]['rescale_GABA_g_intra'] for state in network_state_dict.keys()] # (Baghdoyan, 2012) reference for GABA rescaling

    parameter_space_NMDA_L6_VPM  = [network_state_dict[state]['NMDA_L6_VPM']          for state in network_state_dict.keys()]
    parameter_space_NMDA_L6_TRN  = [network_state_dict[state]['NMDA_L6_TRN']          for state in network_state_dict.keys()]

    for state in network_state_dict.keys():
        if 'w_CT_VPM' in network_state_dict[state].keys():  parameter_space_L6A_VPM     = [network_state_dict[state]['w_CT_VPM']]
        else:                                               parameter_space_L6A_VPM     = [0.75]     # default

        if 'w_CT_TRN' in network_state_dict[state].keys():  parameter_space_L6A_TRN     = [network_state_dict[state]['w_CT_TRN']]
        else:                                               parameter_space_L6A_TRN     = [0.75]     # default
    
    parameter_space_MLe_VPM     = [network_state_dict[state]['w_MLe_VPM']             for state in network_state_dict.keys()]
    parameter_space_VPM_TRN     = [network_state_dict[state]['w_VPM_TRN']             for state in network_state_dict.keys()]
    parameter_space_TRNe_TRNe   = [network_state_dict[state]['gap_gmax']              for state in network_state_dict.keys()]
    
    ##################################################################################################################################################################
    # CT and BS stim datasets
    ct_spk_files_dicts          = [network_state_dict[state]['ct_dict']               for state in network_state_dict.keys()]
    deflection_dataset_path     = [network_state_dict[state]['bs_dict']               for state in network_state_dict.keys()]

    ct_down_times               = [network_state_dict[state]['ct_downTimes']          for state in network_state_dict.keys()]

    ##################################################################################################################################################################

    addHoldingCurrent   =[True]
    addThresholdCurrent =[True]
    ref_rmp_VPM = -58 # mV
    ref_rmp_TRN = -63 # mV
    voltage_step_VPM = 1 # mV
    voltage_step_TRN = 1 # mV
    
    # bWgt_9228_net_16sec  - vecplay      
    parameter_space_shift_targetRMP_VPM = [network_state_dict[state]['shift_RPM'][0]  for state in network_state_dict.keys()]
    parameter_space_shift_targetRMP_TRN = [network_state_dict[state]['shift_RPM'][1]  for state in network_state_dict.keys()]
    group_currents=True

    ############################################################################################################################################################################
    
    target_threshold_currents = []

    calibrated_currents_files=[ '../cells/calibrated_currents/calibrated_currents.json',
                                '../cells/calibrated_currents/calibrated_currents2.json',
                                '../cells/calibrated_currents/calibrated_currents3.json',
                                '../cells/calibrated_currents/calibrated_currents4.json',
                                ]
    map_ThresholdCurrent={}
    for calibrated_currents_file in calibrated_currents_files:
        try:
            with open(calibrated_currents_file, 'r') as file: map_ThresholdCurrent_temp = json.load(file)
            for new_cell in map_ThresholdCurrent_temp.keys():
                if new_cell not in map_ThresholdCurrent.keys():map_ThresholdCurrent.update({new_cell:map_ThresholdCurrent_temp[new_cell]})
        except:
            print('Failed to load file ',calibrated_currents_file) 

    if group_currents:
            
        for ind in range(len(parameter_space_shift_targetRMP_VPM)):
            parameter_space_shift_targetRMP_VPM[ind]

            target_VPL_TC   = ref_rmp_VPM-(parameter_space_shift_targetRMP_VPM[ind]*voltage_step_VPM)
            target_Rt_RC    = ref_rmp_TRN-(parameter_space_shift_targetRMP_TRN[ind]*voltage_step_TRN)
            print('TC: ',target_VPL_TC,' | RT: ',target_Rt_RC)

            map_targetRMP_byCellTemplate = {}
            for cellTemplate in map_ThresholdCurrent.keys():
                if   'VPL_TC' in cellTemplate: map_targetRMP_byCellTemplate.update({cellTemplate:target_VPL_TC})
                elif 'Rt_RC'  in cellTemplate: map_targetRMP_byCellTemplate.update({cellTemplate:target_Rt_RC})

            addThresholdCurrentPops_byCellTemplate={cellTemplate:NetPyNE_BBP.Utils.cubic_extrapolation(map_ThresholdCurrent[cellTemplate]['rmp'], map_ThresholdCurrent[cellTemplate]['i_thresh'], map_targetRMP_byCellTemplate[cellTemplate]) for cellTemplate in map_ThresholdCurrent.keys()}
            target_threshold_currents.append(addThresholdCurrentPops_byCellTemplate)

    else:
        for shift_VPM in parameter_space_shift_targetRMP_VPM: 
            for shift_TRN in parameter_space_shift_targetRMP_TRN: 

                target_VPL_TC   = ref_rmp_VPM-(shift_VPM*voltage_step_VPM)
                target_Rt_RC    = ref_rmp_TRN-(shift_TRN*voltage_step_TRN)

                map_targetRMP_byCellTemplate = {}
                for cellTemplate in map_ThresholdCurrent.keys():
                    if   'VPL_TC' in cellTemplate: map_targetRMP_byCellTemplate.update({cellTemplate:target_VPL_TC})
                    elif 'Rt_RC'  in cellTemplate: map_targetRMP_byCellTemplate.update({cellTemplate:target_Rt_RC})

                addThresholdCurrentPops_byCellTemplate={cellTemplate:NetPyNE_BBP.Utils.cubic_extrapolation(map_ThresholdCurrent[cellTemplate]['rmp'], map_ThresholdCurrent[cellTemplate]['i_thresh'], map_targetRMP_byCellTemplate[cellTemplate]) for cellTemplate in map_ThresholdCurrent.keys()}
                target_threshold_currents.append(addThresholdCurrentPops_byCellTemplate)

    #######################################################################################################################################################################
    #   Vec Play modifications - Not being used
    vp_dict={   'rmp'       :{    'osc':{'VPM':-66, 'TRN':-67},   'awk':{'VPM':-66, 'TRN':-67},   'wxw':{'VPM':-64, 'TRN':-61},     'wxw_weak':{'VPM':-64, 'TRN':-61}},
                'gaba_intra':{    'osc':1.75,                     'awk':0.5,                      'wxw':1.75,                       'wxw_weak':1.0},
                'gaba_extra':{    'osc':1.75,                     'awk':0.1,                      'wxw':1.75,                       'wxw_weak':1.0}}

    # --- GABA weight
    modifyConnWeight_vecPlay_dict ={
                                    'conn|TRN|VPM|chem':                        [(1,     vp_dict['gaba_extra']['awk']),
                                                                                 (4000,  vp_dict['gaba_extra']['wxw_weak']),
                                                                                 (9000,  vp_dict['gaba_extra']['wxw']),
                                                                                 (13000, vp_dict['gaba_extra']['awk']),
                                                                                 ],
                                    'conn|TRN@biophysical|TRN@biophysical|chem':[(1,     vp_dict['gaba_intra']['awk']),
                                                                                 (4000,  vp_dict['gaba_intra']['wxw_weak']),
                                                                                 (9000,  vp_dict['gaba_intra']['wxw']),
                                                                                 (13000, vp_dict['gaba_intra']['awk']),
                                                                                 ], 
                                    }

    #     _wGABA flags
    vector_play_flag_wGABA = [False for key in network_state_dict.keys()]
    for i_key,key in enumerate(network_state_dict.keys()):
        if ('vp_wGABA' in key): vector_play_flag_wGABA[i_key]=True
    #     _wGABA dicts
    vector_play_dict_wGABA = [{} for key in network_state_dict.keys()]
    for flag_wGABA_ind,flag_wGABA in enumerate(vector_play_flag_wGABA): 
        if flag_wGABA: vector_play_dict_wGABA[flag_wGABA_ind]=modifyConnWeight_vecPlay_dict

    # bWgt_9228_net_16sec
    # --- RMP
    shift_t                 = [1, 4000, 9000, 13000]
    target_VPL_TC_shiftRMP  = [vp_dict['rmp']['awk']['VPM'], vp_dict['rmp']['wxw_weak']['VPM'], vp_dict['rmp']['wxw']['VPM'], vp_dict['rmp']['awk']['VPM']]
    target_Rt_RC_shiftRMP   = [vp_dict['rmp']['awk']['TRN'], vp_dict['rmp']['wxw_weak']['TRN'], vp_dict['rmp']['wxw']['TRN'], vp_dict['rmp']['awk']['TRN']]

    modifyRMP_vecPlay_dict={'t':shift_t}    
    modifyRMP_vecPlay_dict.update({key:[] for key in map_ThresholdCurrent.keys()})
    for i in range(len(shift_t)):
        # store_new_targetI={}
        for cellTemplate in map_ThresholdCurrent.keys():
            if   'VPL_TC' in cellTemplate: new_targetRMP = target_VPL_TC_shiftRMP[i]
            elif 'Rt_RC'  in cellTemplate: new_targetRMP = target_Rt_RC_shiftRMP[i]
            modifyRMP_vecPlay_dict[cellTemplate].append(NetPyNE_BBP.Utils.cubic_extrapolation(map_ThresholdCurrent[cellTemplate]['rmp'], map_ThresholdCurrent[cellTemplate]['i_thresh'], new_targetRMP))

    #     RMP flags
    vector_play_flag_RMP = [False for key in network_state_dict.keys()]
    for i_key,key in enumerate(network_state_dict.keys()):
        if ('vp_RMP' in key): vector_play_flag_RMP[i_key]=True
    #     RMP dicts
    vector_play_dict_RMP = [{} for key in network_state_dict.keys()]
    for flag_RMP_ind,flag_RMP in enumerate(vector_play_flag_RMP): 
        if flag_RMP: vector_play_dict_RMP[flag_RMP_ind]=modifyRMP_vecPlay_dict
    
    ############################################################################################################################################################################
    # --- Adding different connectivity schemas between VPM and TRN
    base_dir=os.path.expanduser("~")
    NetPyNE_rootFolder = base_dir+'/Research/Models/BBP/thalamus_netpyne'

    # === Sleep-Wake === #
    # --- barreloid_batch_fig02_sleep_wake
    if      batch_label == 'barreloid_batch_fig02_sleep_wake':
        conn_schema_flags=[                        
                            ['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],     #   Uniform     feedback (convergence)                      (v)
                            ]
    
    # === Awake === #                       and        # === Sleep === #
    # --- barreloid_batch_fig05_awake       and        # --- barreloid_batch_fig06_sleep
    elif    (batch_label == 'barreloid_batch_fig05_awake') or (batch_label == 'barreloid_batch_fig06_sleep'):
        conn_schema_flags=[
                            # --- Open-loop
                            ['VPM_ff_closed' , 'TRN_fb_open' , 'L6A_fb_uniform' ],
                            # --- Mixed Open-loop + Uniform
                            ['VPM_ff_remixed', 'TRN_fb_openRemixed', 'L6A_fb_uniform_tr_75_25'],
                            ['VPM_ff_remixed', 'TRN_fb_openRemixed', 'L6A_fb_uniform_tr_50_50'],
                            ['VPM_ff_remixed', 'TRN_fb_openRemixed', 'L6A_fb_uniform_tr_25_75'],
                            # --- Uniform conns
                            ['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],
                            # --- Mixed Closed-loop + Uniform
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_25_75'],
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_50_50'],
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_75_25'],
                            # --- Closed-loop
                            ['VPM_ff_closed' , 'TRN_fb_closed' , 'L6A_fb_uniform' ],
                            ]
        
    # === Awake CT modulation === #
    # --- barreloid_batch_fig07e_0_uniform_ct_modulation_0   /   barreloid_batch_fig07e_0_uniform_ct_modulation_10   /    barreloid_batch_fig07e_0_uniform_ct_modulation_25   /    barreloid_batch_fig07e_0_uniform_ct_modulation_50   /    barreloid_batch_fig07e_0_uniform_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07e_0_uniform_ct_modulation_100'):
        conn_schema_flags=[
                            ['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],     #   Uniform     feedback (convergence)                      (v)
                            ]
    # --- barreloid_batch_fig07e_1_mixed_ct_modulation_0     /   barreloid_batch_fig07e_1_mixed_ct_modulation_10     /    barreloid_batch_fig07e_1_mixed_ct_modulation_25     /    barreloid_batch_fig07e_1_mixed_ct_modulation_50     /    barreloid_batch_fig07e_1_mixed_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07e_1_mixed_ct_modulation_100'):
        conn_schema_flags=[
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_remixed_tr_50_50'],     #   Mixed-loop  feedback (on-center + convergence)  (new)   (v)
                            ]
    # --- barreloid_batch_fig07e_2_closed_ct_modulation_0    /   barreloid_batch_fig07e_2_closed_ct_modulation_10    /    barreloid_batch_fig07e_2_closed_ct_modulation_25    /    barreloid_batch_fig07e_2_closed_ct_modulation_50    /    barreloid_batch_fig07e_2_closed_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07e_2_closed_ct_modulation_100'):
        conn_schema_flags=[
                            ['VPM_ff_closed' , 'TRN_fb_closed' , 'L6A_fb_closed' ],     #   Closed-loop feedback (on-center)                        (v)
                            ]

    # --- barreloid_batch_fig07f_0_uniform_ct_modulation_0   /   barreloid_batch_fig07f_0_uniform_ct_modulation_10   /    barreloid_batch_fig07f_0_uniform_ct_modulation_25   /    barreloid_batch_fig07f_0_uniform_ct_modulation_50   /    barreloid_batch_fig07f_0_uniform_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07f_0_uniform_ct_modulation_100'):
        conn_schema_flags=[
                            ['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],     #   Uniform     feedback (convergence)                      (v)
                            ]
    # --- barreloid_batch_fig07f_1_mixed_ct_modulation_0     /   barreloid_batch_fig07f_1_mixed_ct_modulation_10     /    barreloid_batch_fig07f_1_mixed_ct_modulation_25     /    barreloid_batch_fig07f_1_mixed_ct_modulation_50     /    barreloid_batch_fig07f_1_mixed_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07f_1_mixed_ct_modulation_100'):
        conn_schema_flags=[
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_50_50'],     #   Mixed-loop  feedback (on-center + convergence)  (new)   (v)
                            ]
    # --- barreloid_batch_fig07f_2_closed_ct_modulation_0    /   barreloid_batch_fig07f_2_closed_ct_modulation_10    /    barreloid_batch_fig07f_2_closed_ct_modulation_25    /    barreloid_batch_fig07f_2_closed_ct_modulation_50    /    barreloid_batch_fig07f_2_closed_ct_modulation_100
    elif    (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_0') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_10') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_25') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_50') or (batch_label == 'barreloid_batch_fig07f_2_closed_ct_modulation_100'):
        conn_schema_flags=[
                            ['VPM_ff_closed' , 'TRN_fb_closed' , 'L6A_fb_uniform' ],     #   Closed-loop feedback (on-center)                        (v)
                            ]

    # === test_TRN_RMP === #
    elif      batch_label == 'test_TRN_RMP':
        # conn_schema_flags=[                        
        #                     ['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],     #   Uniform     feedback (convergence)                      (v)
        #                     ]    
        conn_schema_flags=[
                            # --- Open-loop
                            ['VPM_ff_closed' , 'TRN_fb_open' , 'L6A_fb_uniform' ],
                            # --- Mixed Open-loop + Uniform
                            ['VPM_ff_remixed', 'TRN_fb_openRemixed', 'L6A_fb_uniform_tr_75_25'],
                            ['VPM_ff_remixed', 'TRN_fb_openRemixed', 'L6A_fb_uniform_tr_50_50'],
                            ['VPM_ff_remixed', 'TRN_fb_openRemixed', 'L6A_fb_uniform_tr_25_75'],
                            # --- Uniform conns
                            ['VPM_ff_uniform', 'TRN_fb_uniform', 'L6A_fb_uniform'],
                            # --- Mixed Closed-loop + Uniform
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_25_75'],
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_50_50'],
                            ['VPM_ff_remixed', 'TRN_fb_remixed', 'L6A_fb_uniform_tr_75_25'],
                            # --- Closed-loop
                            ['VPM_ff_closed' , 'TRN_fb_closed' , 'L6A_fb_uniform' ],
                            ]

    _addVersion='_v04'

    conn_schemas = [NetPyNE_rootFolder+'/network_template/barreloid_network_cell_properties__'+schema_ff+'__'+schema_fb+'__'+schema_ct_fb+_addVersion+'.json' for [schema_ff,schema_fb, schema_ct_fb] in conn_schema_flags]
    
    ############################################################################################################################################################################
    
    parameter_space_ih_shift_Thalamus = [0]

    ############################################################################################################################################################################
    
    conversion_factor = 0.001 # must be added to all synapses, except gap junctions

    scale_weights = 1
    scale_EE = conversion_factor * scale_weights * 1.0
    scale_IE = conversion_factor * scale_weights * 1.0
    scale_EI = conversion_factor * scale_weights * 1.0
    scale_II = conversion_factor * scale_weights * 1.0
    
    scale_ML = conversion_factor * scale_weights * 1.0

    weights_dict={
                    'MLe_VPM':      scale_ML * np.array(parameter_space_MLe_VPM),
    
                    'TRNe_TRNe':    np.array(parameter_space_TRNe_TRNe),

                    'TRN_VPM':      scale_IE * np.array(parameter_space_TRN_VPM),
                    'VPM_TRN':      scale_EI * np.array(parameter_space_VPM_TRN),

                    'TRN_TRN':      scale_II * np.array(parameter_space_TRN_TRN),

                    'L6A_VPM':      scale_EE * np.array(parameter_space_L6A_VPM),
                    'L6A_TRN':      scale_EI * np.array(parameter_space_L6A_TRN),
                    
                    # 'L6A_VPM':      1,
                    # 'L6A_TRN':      1,

                    'VPM_L6A':      1,
                    }

    # range variables
    MLe_VPM__weight         = weights_dict['MLe_VPM']
    TRNe_TRNe__weight       = weights_dict['TRNe_TRNe']
    
    TRN_VPM__weight         = weights_dict['TRN_VPM']
    VPM_TRN__weight         = weights_dict['VPM_TRN']
    
    TRN_TRN__weight         = weights_dict['TRN_TRN']
    
    # fixed variables
    L6A_VPM__weight         = weights_dict['L6A_VPM']
    L6A_TRN__weight         = weights_dict['L6A_TRN']
    VPM_L6A__weight         = weights_dict['VPM_L6A']
    
    BaselineDriver_weights_VPM_TRN_ring = [conversion_factor * 1/6] # temporary patch, because weight = 2 was good enough

    ####################################################################################################################################################################################

    # --- Debug sim
    if   session_type == 'debug':       duration = 16000; dt=0.5; sim_min_duration = "0:30:00"; cores = 64; sim_min_duration_expanse = "0:20:00"; cores_expanse = 128
    # --- Full scale sim
    elif session_type == 'detailed':    duration = 16000; dt=0.025; sim_min_duration = "2:35:00"; cores = 64; sim_min_duration_expanse = "2:45:00"; cores_expanse = 64  # bWgt_9401dd4 / bWgt_9401ee1 / bWgt_9401ee2
    # --- (Default) Full scale sim
    else:                               duration = 16000; dt=0.025; sim_min_duration = "2:35:00"; cores = 64; sim_min_duration_expanse = "2:45:00"; cores_expanse = 64  # bWgt_9401dd4 / bWgt_9401ee1 / bWgt_9401ee2
    
    recordStep = 0.1
    if dt>recordStep: recordStep=dt

    recordCurrents = False
    # recordCurrents = True
    if recordCurrents:
        # distributed saving of files during recording of currents
        distributedSaving=True
        recordTraces   = {  "V_soma":                       { "loc": 0.5,   "sec": "soma_0",    "var": "v",},
                            "i__soma_0__SK_E2__ik":         { "loc": 0.5,   "sec": "soma_0",    "var": "ik",    "mech": "SK_E2",},
                            "i__soma_0__TC_HH__ik":         { "loc": 0.5,   "sec": "soma_0",    "var": "ik",    "mech": "TC_HH",},
                            "i__soma_0__TC_HH__ina":        { "loc": 0.5,   "sec": "soma_0",    "var": "ina",   "mech": "TC_HH",},
                            "i__soma_0__TC_Nap_Et2__ina":   { "loc": 0.5,   "sec": "soma_0",    "var": "ina",   "mech": "TC_Nap_Et2",},
                            "i__soma_0__TC_iA__ik":         { "loc": 0.5,   "sec": "soma_0",    "var": "ik",    "mech": "TC_iA",},
                            "i__soma_0__TC_iL__ica":        { "loc": 0.5,   "sec": "soma_0",    "var": "ica",   "mech": "TC_iL",},
                            "i__soma_0__TC_iT_Des98__ica":  { "loc": 0.5,   "sec": "soma_0",    "var": "ica",   "mech": "TC_iT_Des98",},
                            "i__soma_0__TC_ih_Bud97__ih":   { "loc": 0.5,   "sec": "soma_0",    "var": "ih",    "mech": "TC_ih_Bud97",},
                            }
    else:
        # distributed saving of files during recording of currents
        distributedSaving=False
        recordTraces   = {  "V_soma":                       { "loc": 0.5,   "sec": "soma_0",    "var": "v",},}

    initCfg = {
                'duration':         duration,
                'dt':               dt,
                'recordStep':       recordStep,
                'recordTraces':     recordTraces,
                'recordCurrents':   recordCurrents,
                'delayCTvirtual':   0,
                'delayBS':          0,
                ('analysis','plotRaster','timeRange'):[1500,duration],
                ('analysis','plotTraces','timeRange'):[1500,duration],
                ('analysis','plotLFP',   'timeRange'):[1500,duration],
                'distributedSaving':distributedSaving,
                }
    if (batch_label == 'barreloid_batch_fig06_sleep'): initCfg.update({'record_windowed_spectrogram':True})
    
    params = {  
                # --- Grouped params
                'deflection_dataset_path':                  deflection_dataset_path,
                'addThresholdCurrentPops_byCellTemplate':   target_threshold_currents,
                'ct_spk_files_dict':                        ct_spk_files_dicts,
                'ct_downTimes':                             ct_down_times,

                'replace_spk_file':                         parameter_space_replace_spk_file, # (2024_11_10) - replacing spikes in CT virtual by artificial angular tuned spikes

                'add_KLeak__VPM':                           parameter_space_KLeak_VPM,
                'add_KLeak__TRN':                           parameter_space_KLeak_TRN,
                'TRN_VPM__weight':                          TRN_VPM__weight,
                'TRN_TRN__weight':                          TRN_TRN__weight,
                'TRNe_TRNe__weight':                        TRNe_TRNe__weight,
                'VPM_TRN__weight':                          VPM_TRN__weight,
                'MLe_VPM__weight':                          MLe_VPM__weight,

                'modifyNMDAratio_L6A_VPM':                  parameter_space_NMDA_L6_VPM,
                'modifyNMDAratio_L6A_TRN':                  parameter_space_NMDA_L6_TRN,

                #  -  cell params to improve neuromodulation
                'rescale_iNap': parameter_space_rescale_inap_gmax,     # (default 1)
                'rescale_ih':   parameter_space_rescale_ih_gmax,       # (default 1)
                'rescale_iT':   parameter_space_rescale_it_gmax,       # (default 1)
                'rescale_iA':   parameter_space_rescale_ia_gmax,       # (default 1)

                #  -  params to add vecplay modification
                'addModifyRMP_vecPlay':                     vector_play_flag_RMP,
                'addModifyConnWeight_vecPlay':              vector_play_flag_wGABA,
                'modifyRMP_vecPlay_dict':                   vector_play_dict_RMP,
                'modifyConnWeight_vecPlay_dict':            vector_play_dict_wGABA,
                
                # --- Range params
                'dump_cell_properties':                     conn_schemas,
                'add_iT_shift__Thalamus':                   parameter_space_ih_shift_Thalamus, # shifts iT for all thalamic cells
                'addUSE_MLe':                               parameter_space_addUSE_MLe,
                'rescaleDep_MLe':                           parameter_space_rescaleDep_MLe,
                
                # --- Single value params
                'L6A_VPM__weight':                          L6A_VPM__weight,
                'L6A_TRN__weight':                          L6A_TRN__weight,
                'addHoldingCurrent':                        addHoldingCurrent,
                'addThresholdCurrent':                      addThresholdCurrent,
                'BaselineDriver_weights_VPM_TRN_ring':      BaselineDriver_weights_VPM_TRN_ring,

                'rescale_ih_gmax__VPM':                     parameter_space_rescale_iH_gmax_VPM,
                'rescale_ih_gmax__TRN':                     parameter_space_rescale_iH_gmax_TRN,
                }

    b = Batch(cfgFile='barreloid_cfg.py', netParamsFile='barreloid_netParams.py', params=params, initCfg=initCfg,
              groupedParams=[                                
                                'deflection_dataset_path',
                                'addThresholdCurrentPops_byCellTemplate',
                                'ct_spk_files_dict',
                                'ct_downTimes',

                                'replace_spk_file', # (2024_11_10) - replacing spikes in CT virtual by artificial angular tuned spikes

                                'add_KLeak__VPM',
                                'add_KLeak__TRN',
                                'TRN_VPM__weight',
                                'TRN_TRN__weight',
                                'TRNe_TRNe__weight',
                                
                                'VPM_TRN__weight',
                                'MLe_VPM__weight',

                                'modifyNMDAratio_L6A_VPM',
                                'modifyNMDAratio_L6A_TRN',

                                'addModifyRMP_vecPlay',
                                'addModifyConnWeight_vecPlay',
                                'modifyRMP_vecPlay_dict',
                                'modifyConnWeight_vecPlay_dict',

                                #  -  cell params to improve neuromodulation
                                'rescale_iNap',
                                'rescale_ih',
                                'rescale_iT',
                                'rescale_iA',

                             ])
    
    ########################################################
    
    b.batchLabel = batch_flag+'_'+batch_label

    ########################################################

    b.saveFolder = '../data/batch_sims/'+b.batchLabel
    b.method = 'grid'

    base_dir=os.path.expanduser("~")
    if base_dir == '/home/jmoreira':
        b.runCfg = {'type': 'hpc_slurm',
                    # 'allocation': 'TG-IBN140002',
                    'allocation': 'TG-MED240050',
                    'partition': 'compute', # this partition was throwing submission errors on 2024_08_07 - switched to partition large-shared
                    # 'partition': 'shared',
                    # 'partition': 'large-shared',
                    'walltime': sim_min_duration_expanse,
                    'nodes': 1,
                    'coresPerNode': cores_expanse,
                    'email': 'joao.moreira@downstate.edu',
                    'folder': '/home/jmoreira/Research/Models/BBP/thalamus_netpyne/sim/',
                    'script': 'barreloid_init.py',
                    'mpiCommand': 'mpiexec',
                    'custom': '\n#SBATCH --export=ALL\n#SBATCH --partition=compute',
                    'skipCfg': True,
                    }
    else:
        # Set output folder, grid method (all param combinations), and run configuration
        b.runCfg = {'type': 'hpc_sge',
                    'jobName': b.batchLabel,
                    'cores': cores,
                    'script': 'barreloid_init.py',
                    'walltime': sim_min_duration,
                    'vmem': '128G',
                    'skipCfg': True,
                    }

    b.run()

###################################################################################################################################################################
# Main code
if __name__ == '__main__':

    ##############################################################################################################################################################
    # --- Type of simulation (changing dt to speed up runtime during debug)
    # sessionType = 'debug'
    sessionType = 'detailed'

    # --- Naming flag added before the batchLabel when saving the files
    batchFlag = 'paper_002'

    ##############################################################################################################################################################
    # RUN SIMULATIONS FOR FIGURE 02
    # # === Sleep-Wake === #
    batchLabel = 'barreloid_batch_fig02_sleep_wake' ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    # RUN SIMULATIONS FOR FIGURE 05
    # === Awake === #
    batchLabel = 'barreloid_batch_fig05_awake'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    # RUN SIMULATIONS FOR FIGURE 06
    # === Sleep === #
    batchLabel = 'barreloid_batch_fig06_sleep'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    ## RUN SIMULATIONS FOR FIGURE 07
    # === Awake CT modulation === # - Closed-loop CT
    #  -  Uniform  -  #
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_0'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_10'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_25'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_50'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_0_uniform_ct_modulation_100'    ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Mixed  -  #
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_0'        ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_10'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_25'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_50'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_1_mixed_ct_modulation_100'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Closed  -  #
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_0'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_10'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_25'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_50'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07e_2_closed_ct_modulation_100'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    
    # === Awake CT modulation === # - UNIFORM CT
    #  -  Uniform  -  #
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_0'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_10'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_25'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_50'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_0_uniform_ct_modulation_100'    ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Mixed  -  #
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_0'        ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_10'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_25'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_50'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_1_mixed_ct_modulation_100'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    #  -  Closed  -  #
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_0'       ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_10'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_25'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_50'      ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)
    batchLabel = 'barreloid_batch_fig07f_2_closed_ct_modulation_100'     ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)    
    

    ###################################################################################################################################################################
    # # # # === test_TRN_RMP === #
    # # batchLabel = 'test_TRN_RMP' ;           batchWeight(batch_label = batchLabel, session_type=sessionType, batch_flag=batchFlag)

    # # batchWeight(batch_label = batchLabel)


###################################################################################################################################################################