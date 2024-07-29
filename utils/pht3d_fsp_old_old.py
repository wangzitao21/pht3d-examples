#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is the Python tool "PHT3D-FSP" (PHT3D - FloPy Support Package) v0.12, see readme file for further information.

This work was realized in the Deutsche Forschungsgemeinschaft (DFG) project Reactive Transport (GR 4514/3-1) within the research unit FOR 5094: The dynamic deep subsurface of high-energy beaches (DynaDeep).

This work is published under the GNU GENERAL PUBLIC LICENSE Version 3.

@author: Stephan L. Seibert (stephan.seibert@uol.de) and Janek Greskowiak (janek.greskowiak@uol.de)
"""

import numpy as np
import pandas as pd

def create(xlsx_path="./", xlsx_name="pht3d_species.xlsx",
           nlay=1, nrow=1, ncol=1,
           ph_os=2, ph_temp=25.0, ph_asbin=0,
           ph_eps_aqu=1e-10, ph_ph=1e-3,
           ph_print=0,
           ph_cb_offset=0, ph_surf_calc_type="-diffuse_layer",
           write_ph="yes"
          ):

    """
    Input:
        xlsx_path: Excel 文件路径
        xlsx_name: Excel 文件名
        nlay: 层数
        nrow: 行数
        ncol: 列数
        ph_os: 分裂算子方案
        ph_temp: 反应温度 (摄氏度)。
        ph_asbin: 决定是否生成包含每个网格单元在每个输出时间点的浓度值的 ASCII 输出文件
        ph_eps_aqu: 用来决定是否需要对网格单元进行化学反应计算的阈值, 如果某个网格单元中所有移动物种的浓度变化量都小于该值, 则认为该网格单元的化学反应可以忽略, 从而减少计算量
        ph_ph: 同上阈值, 如果某个网格单元中 pH 的变化量小于该值, 则该网格单元的化学反应可以忽略, 从而减少计算量
        ph_print: 是否打印附加输出, 没懂
        ph_cb_offset: 用来决定是否需要将溶液的电荷不平衡量进行迁移
        ph_surf_calc_type: 表面复合计算的类型
        write_ph: 是否写入到 pht3d_ph.dat 文件
    Output:
        spec_array_dict: 包含所有物种初始浓度数组的字典

        pht3d_fsp.create.sconc_btn: 用于 Flopy BTN 包的初始浓度变量字符串, 输出示例: "sconc = spec['Species'], sconc2 = spec['ph'], sconc3 = spec['pe']"
        pht3d_fsp.create.crch_ssm: 用于 Flopy SSM 包的补给浓度变量字符串, 输出示例: "crch = rch_spec['Species'], crch2 = rch_spec['ph'], crch3 = rch_spec['pe']"
        pht3d_fsp.create.ncomp: 用于模拟的溶质总数
        pht3d_fsp.create.mcomp: 可移动的溶质数量
    """

    # Define all species here...
    species = pd.read_excel(xlsx_path+xlsx_name)
    species_array = np.array(species)    
    
    # Create ditionary with species number and species initial concentration
    spec_dict = {}
    spec_array_dict = {}
    for i in np.arange(0,len(species)):
        key="create." + species.name[i] 
        spec_dict[key] = species.initial_concentration[i]
        spec_array_key = species.name[i]
        spec_array_dict[spec_array_key] = species.initial_concentration[i]*np.ones((nlay, nrow, ncol), dtype=np.float32)
           
    # Create initial concentration arrays for all species
    for key,value in spec_dict.items():
        exec(f'{key} = np.ones((nlay, nrow, ncol), dtype=np.float32)*{value}')
            
    # Create 'sconc' variables for flopy BTN package
    create.sconc_btn = []
    for i in np.arange(0,len(species)):
        if i == 0:
             word1 = str("sconc = " + "spec['" + species.name[i] + "'], ")
        else:
             word1 = str("sconc"+"{:0>1d}".format(i+1)+" = " + "spec['" + species.name[i] + "'], ")  # changed by JG 7.1.23
        create.sconc_btn.append(word1)
    create.sconc_btn = ''.join(create.sconc_btn)
    create.sconc_btn = create.sconc_btn[:-2]
    
    # Create 'crch' variables for flopy SSM package
    create.crch_ssm = []
    for i in np.arange(0,len(species)):
        if i == 0:
            word1 = "crch = " + "rch_spec['" + species.name[i] + "'], " # changed by JG 7.1.23
        else:
            word1 = "crch"+"{:0>1d}".format(i+1) + " = rch_spec['" + species.name[i] + "'], " # changed by JG 7.1.23
        create.crch_ssm.append(word1)
    create.crch_ssm = ''.join(create.crch_ssm)
    create.crch_ssm = create.crch_ssm[:-2]
    
    # Calculate amount of the different reactant types
    mobile_comp = np.nansum(species.mobility == 'mobile')
    immobile_comp = np.nansum(species.mobility == 'immobile')
    create.ncomp = mobile_comp + immobile_comp
    create.mcomp = mobile_comp
    kinetic_mobile_species = np.nansum(species.type == 'A')
    LEA_aqueous_species = np.nansum(species.type == 'B')
    kinetic_immobile_species = np.nansum(species.type == 'C')
    LEA_minerals = np.nansum(species.type == 'D')
    exchange_species = np.nansum(species.type == 'E')
    ion_exchange_no = np.nansum(species.ion_exchange == 'yes')
    surface_species = np.nansum(species.type == 'F')    
    kinetic_minerals = np.nansum(species.type == 'G')
    kinetic_surface_complexation = np.nansum(species.type == 'H')    
    
    # Create lists for the different reactant types
    # 1. Mobile kinetic species (Type 'A')
    kinetic_mobile_species_names = [] # containing name
    kinetic_mobile_species_number = [] # containing total number of kinetic parameters for each kinetic component
    kinetic_mobile_species_parameters = [] # containing the kinetic parameters for each kinetic component
    kinetic_mobile_species_formula = [] # containing the formulas for each kinetic component
    kinetic_mobile_species_m0 = [] # containing the initial amount of component
    for i in np.arange(0,len(species)):
        if species.type[i] == 'A':
            kinetic_mobile_species_names.append(species.species[i])
            kinetic_parameters = []
            kinetic_parameters = species_array[i,14:115]
            kinetic_parameters = kinetic_parameters.astype("float")
            kinetic_parameters = kinetic_parameters[~np.isnan(kinetic_parameters)]
            kinetic_mobile_species_parameters.append(kinetic_parameters)
            kinetic_mobile_species_number.append(len(kinetic_parameters))
            kinetic_mobile_species_formula.append(species.formula[i])
            kinetic_mobile_species_m0.append(float(species.argument[i]))    
    # 2. LEA solution master species (Type 'B')
    LEA_species_names = [] # containing name of kinetic mobile species
    LEA_species_argument = [] # containing the formulas for each kinetic mobile species
    for i in np.arange(0,len(species)):
        if species.type[i] == 'B':
            LEA_species_names.append(species.species[i])
            LEA_species_argument.append(species.argument[i])
    # 3. Immobile kinetic species (Type 'C')
    kinetic_immobile_species_names = [] # containing name
    kinetic_immobile_species_number = [] # containing total number of kinetic parameters for each kinetic component
    kinetic_immobile_species_parameters = [] # containing the kinetic parameters for each kinetic component
    kinetic_immobile_species_formula = [] # containing the formulas for each kinetic component
    kinetic_immobile_species_m0 = [] # containing the initial amount of component
    for i in np.arange(0,len(species)):
        if species.type[i] == 'C':                    
            kinetic_immobile_species_names.append(species.species[i])
            kinetic_parameters = []
            kinetic_parameters = species_array[i,14:115]
            kinetic_parameters = kinetic_parameters.astype("float")
            kinetic_parameters = kinetic_parameters[~np.isnan(kinetic_parameters)]
            kinetic_immobile_species_parameters.append(kinetic_parameters)
            kinetic_immobile_species_number.append(len(kinetic_parameters))
            kinetic_immobile_species_formula.append(species.formula[i])
            kinetic_immobile_species_m0.append(float(species.argument[i]))
    # 4. LEA minerals (Type 'D')
    LEA_minerals_names = [] # containing name of kinetic mobile species
    LEA_minerals_argument = [] # containing the formulas for each kinetic mobile species
    for i in np.arange(0,len(species)):
        if species.type[i] == 'D':
            LEA_minerals_names.append(species.species[i])
            LEA_minerals_argument.append(species.argument[i])
    # 5. Exchange species (Type 'E')
    ion_exchange_species_names = [] # containing name of kinetic mobile species
    ion_exchange_species_stoichiometry = [] # containing stoichiometry of species wrt exchange master species
    ion_exchange_species_master_species = [] # containing the exchange master species
    ion_exchange_species_type = [] # containing the species type
    for i in np.arange(0,len(species)):
        if species.ion_exchange[i] == 'yes':
            ion_exchange_species_names.append(species.species[i])
            ion_exchange_species_stoichiometry.append(species.exchange_stoichiometry[i])
            ion_exchange_species_master_species.append(species.exchange_master_species[i])
            ion_exchange_species_type.append(species.type[i])
    # 6. Surface master species (Type 'F')
    surface_complexation_names = [] # containing name of surface binding site
    surface_complexation_area = [] # containing specific surface area
    surface_complexation_mass = [] # containing mass of solid to calculate surface area
    surface_complexation_phase = [] # containing pure phase or kinetic reactant, copuled to the binding site
    surface_complexation_switch = [] # containing argument declaring type of phase
    for i in np.arange(0,len(species)):
        if species.type[i] == 'F':
            surface_complexation_names.append(species.species[i])
            surface_complexation_area.append(species.surface_area[i])
            surface_complexation_mass.append(species.surface_mass[i])
            surface_complexation_phase.append(species.surface_phase[i])
            surface_complexation_switch.append(species.surface_switch[i])
    # 7. Kinetic minerals (Type 'G')
    kinetic_minerals_names = [] # containing name
    kinetic_minerals_number = [] # containing total number of kinetic parameters for each kinetic component
    kinetic_minerals_parameters = [] # containing the kinetic parameters for each kinetic component
    kinetic_minerals_formula = [] # containing the formulas for each kinetic component
    for i in np.arange(0,len(species)):
        if species.type[i] == 'G':
            kinetic_minerals_names.append(species.species[i])
            kinetic_parameters = []
            kinetic_parameters = species_array[i,14:115]
            kinetic_parameters = kinetic_parameters.astype("float")
            kinetic_parameters = kinetic_parameters[~np.isnan(kinetic_parameters)]
            kinetic_minerals_parameters.append(kinetic_parameters)
            kinetic_minerals_number.append(len(kinetic_parameters))
            kinetic_minerals_formula.append(species.formula[i])
    # 8. Kinetic surface complexation (Type 'H')
    kinetic_suface_complexation_names = [] # containing name
    kinetic_suface_complexation_number = [] # containing total number of kinetic parameters for each kinetic component
    kinetic_suface_complexation_parameters = [] # containing the kinetic parameters for each kinetic component
    kinetic_suface_complexation_formula = [] # containing the formulas for each kinetic component
    for i in np.arange(0,len(species)):
        if species.type[i] == 'H':
            kinetic_suface_complexation_names.append(species.species[i])
            kinetic_parameters = []
            kinetic_parameters = species_array[i,14:115]
            kinetic_parameters = kinetic_parameters.astype("float")
            kinetic_parameters = kinetic_parameters[~np.isnan(kinetic_parameters)]
            kinetic_suface_complexation_parameters.append(kinetic_parameters)
            kinetic_suface_complexation_number.append(len(kinetic_parameters))
            kinetic_suface_complexation_formula.append(species.formula[i])      
    
    #%% Create 'pht3d_ph.dat'
    
    ### Define global parameters
    # PH1 Record
    l01_01_OS = ph_os # operator-splitting scheme
    l01_02_TEMP = ph_temp # temperature [°C] for reactions with temperature dependence; negative number refers to species representing temperature
    l01_03_ASBIN = ph_asbin # format of concentration outputs for all cells and stress periods (0=binary,1=binary&ASCII)
    l01_04_EPS_AQU = ph_eps_aqu# minimum amount [moles] of change for a species between reaction steps required for calls to PHREEQC-2
    l01_05_EPS_PH = ph_ph # minimum change of pH between reaction steps required for calls to PHREEQC-2
    l01_06_PRINT = ph_print # print additional PHREEQC output; not described in PHT3D manual
    
    # PH2 Record
    l02_01_CB_OFFSET = ph_cb_offset # flag to indicate if charge imbalance of solution is to be transported; default 0
    
    # PH3 Record
    l03_01_NR_SOL_MST_SPEC_EQU = LEA_aqueous_species # number of aqueous LEA species (incl. pH and pe)
    
    # PH4 Record
    l04_01_NR_MIN_EQU = LEA_minerals # number of (immobile) LEA minerals
    
    # PH5 Record
    l05_01_NR_ION_EX = exchange_species # number of exchange species
    
    # PH6 Record
    l06_01_NR_SURF = surface_species # number of surface master species
    
    # PH7 Record
    l07_01_NR_MOB_KIN = kinetic_mobile_species # number of kinetic, non-LEA mobile species with rate expression
    l07_02_NR_MIN_KIN = kinetic_minerals # number of kinetic, non-LEA minerals with rate expression
    l07_03_NR_SURF_KIN = kinetic_surface_complexation # number of kinetic surface-complexation reactions; not supported by 'pht3d_ph.dat' yet, thus, must be set to '0'
    l07_04_NR_SURF_KIN = kinetic_immobile_species # number of kinetic, non-LEA immobile species with rate expression
    
    # PH13 Record
    SURF_CALC_TYPE = ph_surf_calc_type # define type of SCM calculation used for surface complexation (provide empty string if no option shall be used)
    
    ### Write 'pht3d_ph.dat'
    if write_ph == "yes":
        with open("./pht3d_ph.dat",'w') as pht3d_ph:
            # PH1 Record:
            pht3d_ph.write(str(l01_01_OS)+' '+str(l01_02_TEMP)+' '+str(l01_03_ASBIN)+' '+str(l01_04_EPS_AQU)+' '+str(l01_05_EPS_PH)+' '+str(l01_06_PRINT)+'\n')
            # PH2 Record:
            pht3d_ph.write(str(l02_01_CB_OFFSET)+'\n')
            # PH3 Record:
            pht3d_ph.write(str(l03_01_NR_SOL_MST_SPEC_EQU)+'\n')
            # PH4 Record:
            pht3d_ph.write(str(l04_01_NR_MIN_EQU)+'\n')
            # PH5 Record:
            pht3d_ph.write(str(l05_01_NR_ION_EX)+'\n')
            # PH6 Record:
            pht3d_ph.write(str(l06_01_NR_SURF)+'\n')
            # PH7 Record:
            pht3d_ph.write(str(l07_01_NR_MOB_KIN)+' '+str(l07_02_NR_MIN_KIN)+' '+str(l07_03_NR_SURF_KIN)+' '+str(l07_04_NR_SURF_KIN)+'\n')
            # PH8 Record:
            if kinetic_mobile_species > 0:
                for i in np.arange(0,kinetic_mobile_species):
                    argument_value = []
                    argument_value = str(kinetic_mobile_species_m0[i])
                    argument_value = argument_value.lower()
                    if argument_value == 'nan':
                        pht3d_ph.write(str(kinetic_mobile_species_names[i])+' '+str(kinetic_mobile_species_number[i])+'\n')
                    else:
                        pht3d_ph.write(str(kinetic_mobile_species_names[i])+' '+str(kinetic_mobile_species_number[i])+' '+str(kinetic_mobile_species_m0[i])+'\n')
                    if kinetic_mobile_species_number[i] > 0:
                        for j in np.arange(0,len(kinetic_mobile_species_parameters[i])):
                            pht3d_ph.write(str(kinetic_mobile_species_parameters[i][j])+'\n')
                    pht3d_ph.write('-Formula '+str(kinetic_mobile_species_formula[i])+'\n')
            # PH9 Record:
            for i in np.arange(0,LEA_aqueous_species):
                argument_value = []
                argument_value = str(LEA_species_argument[i])
                argument_value = argument_value.lower()
                if argument_value == 'nan':
                    pht3d_ph.write(str(LEA_species_names[i])+'\n')
                else:
                    pht3d_ph.write(str(LEA_species_names[i])+' '+str(LEA_species_argument[i])+'\n')
            # PH10 Record:
            if kinetic_immobile_species > 0:
                for i in np.arange(0,kinetic_immobile_species):
                    argument_value = []
                    argument_value = str(kinetic_immobile_species_m0[i])
                    argument_value = argument_value.lower()
                    if argument_value == 'nan':
                        pht3d_ph.write(str(kinetic_immobile_species_names[i])+' '+str(kinetic_immobile_species_number[i])+'\n')
                    else:
                        pht3d_ph.write(str(kinetic_immobile_species_names[i])+' '+str(kinetic_immobile_species_number[i])+' '+str(kinetic_immobile_species_m0[i])+'\n')
                    if kinetic_immobile_species_number[i] > 0:
                        for j in np.arange(0,len(kinetic_immobile_species_parameters[i])):
                            pht3d_ph.write(str(kinetic_immobile_species_parameters[i][j])+'\n')
                    pht3d_ph.write('-Formula '+str(kinetic_immobile_species_formula[i])+'\n')
            # PH11 Record:
            if LEA_minerals > 0:
                for i in np.arange(0,LEA_minerals):
                    argument_value = []
                    argument_value = str(LEA_minerals_argument[i])
                    argument_value = argument_value.lower()
                    if argument_value == 'nan':
                        pht3d_ph.write(str(LEA_minerals_names[i])+'\n')
                    elif argument_value != 'nan':
                        pht3d_ph.write(str(LEA_minerals_names[i])+' '+str(LEA_minerals_argument[i])+'\n')
            # PH12 Record:
            if exchange_species > 0:
                for i in np.arange(0,len(ion_exchange_species_names)):
                    argument_value = []
                    argument_value = str(ion_exchange_species_master_species[i])
                    argument_value = argument_value.lower()
                    if ion_exchange_species_type[i] == 'E':
                        pht3d_ph.write(str(ion_exchange_species_names[i])+'\n')
                    elif argument_value != 'nan':
                        pht3d_ph.write(str(ion_exchange_species_names[i])+' '+str(int(ion_exchange_species_stoichiometry[i]))+' '+str(ion_exchange_species_master_species[i])+'\n')
                    else:
                        pht3d_ph.write(str(ion_exchange_species_names[i])+' '+str(int(ion_exchange_species_stoichiometry[i]))+'\n')
            # PH13 Record:
            if surface_species > 0:
                for i in np.arange(0,len(surface_complexation_names)):
                    argument_value = []
                    argument_value = str(surface_complexation_phase[i])
                    argument_value = argument_value.lower()
                    if argument_value == 'nan':
                        pht3d_ph.write(str(surface_complexation_names[i])+' '+str(surface_complexation_area[i])+' '+str(surface_complexation_mass[i])+'\n')
                    else:
                        pht3d_ph.write(str(surface_complexation_names[i])+' '+str(surface_complexation_area[i])+' '+str(surface_complexation_mass[i])+' '+str(surface_complexation_phase[i])+' '+str(surface_complexation_switch[i])+'\n')
                if len(SURF_CALC_TYPE) > 0:
                    pht3d_ph.write(SURF_CALC_TYPE+'\n')
            # PH14 Record:
            if kinetic_minerals > 0:
                for i in np.arange(0,kinetic_minerals):
                    pht3d_ph.write(str(kinetic_minerals_names[i])+' '+str(kinetic_minerals_number[i])+'\n')
                    if kinetic_minerals_number[i] > 0:
                        for j in np.arange(0,len(kinetic_minerals_parameters[i])):
                            pht3d_ph.write(str(kinetic_minerals_parameters[i][j])+'\n')  
                    argument_value = str(kinetic_minerals_formula[i])    
                    argument_value = argument_value.lower()  
                    if argument_value == 'nan':        
                       pht3d_ph.write('-Formula \n')  
                    else:
                       pht3d_ph.write('-Formula '+str(kinetic_minerals_formula[i])+'\n') 
    return spec_array_dict