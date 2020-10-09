#ifndef _H_GENERATION_SPECIFICATIONS_
#define _H_GENERATION_SPECIFICATIONS_

#include <LeptonWeighter/Generator.h>
#include <LeptonWeighter/Event.h>
#include <LeptonWeighter/ParticleType.h>
#include "GolemTools.h"
#include "GolemEnumDefinitions.h"
#include "GolemMCSet.h"

namespace golemfit {
  namespace sterile {
    std::map<std::string,MCSet> GetSimInfo(std::string mc_dataPath)
    {
       std::map<std::string,MCSet> simInfo = {
        /*
        * OFFICIAL ARES GOLDEN SETS
        */
        {"ares_golden_macho_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_98765.h5",{false,0},	98765/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_superheavy_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_50000.h5",{false,0},50000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_heavy_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_25000.h5",{false,0},	25000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_feather_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_10000.h5",{false,0},	10000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_welter_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_5000.h5",{false,0},	5000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_lite_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_500.h5",{false,0},500/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_wimp_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_5.h5",{false,0},5/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_golden_15000.h5",{false,0},15000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES DIAMOND SETS
        */
        {"ares_diamond_macho_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_99630.h5",{false,0},	99630/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_superheavy_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_50000.h5",{false,0},50000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_heavy_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_25000.h5",{false,0},	25000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_feather_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_10000.h5",{false,0},	10000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_welter_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_5000.h5",{false,0},	5000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_lite_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_500.h5",{false,0},	500/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_wimp_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_5.h5",{false,0},	5/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_diamond_15000.h5",{false,0},		15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /* OFFICIAL ARES DIAMOND SPLIT SETS
        */
        {"ares_diamond_welter_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_5000",{true,5},5000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES Platinum
        */
        {"ares_HQE132_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HQE132_HE_platinum_69920.h5",{false,0},	69920/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HQE132_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HQE132_LE_platinum_5755.h5",{false,0},		5755/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HQE138_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HQE138_HE_platinum_60430.h5",{false,0},	60430/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HQE138_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HQE138_LE_platinum_7065.h5",{false,0},		7065/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },


        {"ares_platinum_feather_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_platinum_10000.h5",{false,0},	10000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"BI1_platinum_feather_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"BI1_Nominal_HE_platinum_10000.h5",{false,0},		10000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_welter_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_platinum_5000.h5",{false,0},		5000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_titanium_welter_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_titanium_5000.h5",{false,0},		5000/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_lite_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_platinum_500.h5",{false,0},		500/5,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},		15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_feather_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},		15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_titanium_welter_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_titanium_14915.h5",{false,0},		14915/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /* OFFICIAL ARES PLATINUM SPLIT SETS
        */
        {"ares_platinum_welter_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_5000",{true,5},5000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_feather_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_10000",{true,10},10000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_heavy_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_25000",{true,25},25000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_super_heavy_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_50000",{true,50},50000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_macho_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_95000",{true,95},95000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_mighty_split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_75000",{true,75},75000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_mighty_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_IC86.AVG_1.27_lite_platinum_75000.h5",{false,0},75000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"Ares_Nominal_Skip_HE_platinum_10000",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_Skip_HE_platinum_10000",{true,10},10000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_platinum_welter_split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_feather_split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_heavy_split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_super_heavy_split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_macho_split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_mighty_split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_platinum_mighty_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},15000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"Platinum_Split",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_98000",{true,98},98000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true)))
        },
        {"Platinum_lite",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_lite_500.h5",{false,0},500/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true)))
        },
        {"Platinum_Central_Split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Noise_LE_platinum_14615.h5",{false,0},14615/5.,1.27,-1,1.0,
                    LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"Platinum_Central_Split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Noise_HE_platinum_90000",{true,90},90000/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"Platinum_SpiceMie_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_SpiceMie_LE_18090.h5",{false,0},18090/5.,1.27,-1,1.0,
                    LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"Platinum_SpiceMie_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_SpiceMie_HE_73155.h5",{false,0},73155/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },

        {"Platinum_Central_Feather_Split_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_Central_Feather_Split_LE.h5",{false,0},14655/5.,1.27,-1,1.0,
                    LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"Platinum_Central_Feather_Split_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Noise_HE_platinum_2500.h5",{false,0},2500/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"Platinum_Welter_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_Welter_LE.h5",{false,0},14655/5.,1.27,-1,1.0,
                    LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"Platinum_Welter_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Platinum_Welter_HE.h5",{false,0},2500/5.,1.27,-1,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES DOMEFF SETS
        */
        {"ares_DE_1.23_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.23_HE_diamond_99645.h5",{false,0},	99645/5.,1.23,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.25_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.25_HE_diamond_99665.h5",{false,0},	99665/5.,1.25,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.27_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.27_HE_diamond_99630.h5",{false,0},	99630/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.30_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.30_HE_diamond_99630.h5",{false,0},	99630/5.,1.30,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.33_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.33_HE_diamond_99665.h5",{false,0},	99665/5.,1.33,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.35_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.35_HE_diamond_99510.h5",{false,0},	99510/5.,1.35,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },


        {"ares_DE_1.23_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.23_LE_diamond_15000.h5",{false,0},	15000/5.,1.23,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.25_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.25_LE_diamond_14985.h5",{false,0}, 14985/5.,1.25,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.27_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.27_LE_diamond_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.30_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.30_LE_diamond_15000.h5",{false,0},	15000/5.,1.30,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.33_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.33_LE_diamond_14995.h5",{false,0},	14995/5.,1.33,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.35_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.35_LE_diamond_14975.h5",{false,0},	14975/5.,1.35,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES PLATINUM DOMEFF SETS
        */
        {"ares_DE_1.23_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.23_HE_platinum_99645.h5",{false,0},	99645/5.,1.23,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.25_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.25_HE_platinum_99665.h5",{false,0},	99665/5.,1.25,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.27_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_platinum_99630.h5",{false,0},	99630/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.30_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.30_HE_platinum_99605.h5",{false,0},	99605/5.,1.30,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.33_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.33_HE_platinum_99665.h5",{false,0},	99665/5.,1.33,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.35_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.35_HE_platinum_99510.h5",{false,0},	99510/5.,1.35,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },

        {"ares_DE_1.23_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.23_LE_platinum_15000.h5",{false,0},	15000/5.,1.23,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.25_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.25_LE_platinum_14985.h5",{false,0},	14985/5.,1.25,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.27_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.27_LE_platinum_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.30_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.30_LE_platinum_15000.h5",{false,0},	15000/5.,1.30,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.33_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.33_LE_platinum_14995.h5",{false,0},	14995/5.,1.33,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.35_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.35_LE_platinum_14975.h5",{false,0},	14975/5.,1.35,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES GOLDEN DOMEFF SETS
        */
        {"ares_DE_1.23_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.23_HE_golden_99645.h5",{false,0},	99645/5.,1.23,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.25_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.25_HE_golden_99665.h5",{false,0},	99665/5.,1.25,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.27_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.27_HE_golden_98765.h5",{false,0},	98765/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.30_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.30_HE_golden_99630.h5",{false,0},	99630/5.,1.30,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.33_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.33_HE_golden_99665.h5",{false,0},	99665/5.,1.33,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_DE_1.35_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.35_HE_golden_99510.h5",{false,0},	99510/5.,1.35,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },

        {"ares_DE_1.23_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.23_LE_golden_15000.h5",{false,0},	15000/5.,1.23,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.25_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.25_LE_golden_14985.h5",{false,0},	14985/5.,1.25,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.27_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.27_LE_golden_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.30_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.30_LE_golden_15000.h5",{false,0},	15000/5.,1.30,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.33_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.33_LE_golden_14995.h5",{false,0},	14995/5.,1.33,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_DE_1.35_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_DE_1.35_LE_golden_14975.h5",{false,0},	14975/5.,1.35,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL PLATINUM ARES HOLEICE SETS
        * HI1 = 2.0
        * HI2 = -1.0
        * HI3 = -5.0
        * HI4 = 0.5
        * HI5 = -3.0
        */ 
        {"ares_HI_2.0_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_HE_platinum_99510.h5",{false,0},	99510/5.,1.27,2.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-1.0_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_platinum_99630.h5",{false,0},	99630/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-5.0_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI3_HE_platinum_97840.h5",{false,0},	97840/5.,1.27,-5.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_0.5_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_HE_platinum_95655.h5",{false,0},	95655/5.,1.27,0.5,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-3.0_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI5_HE_platinum_95225.h5",{false,0},	95225/5.,1.27,-3.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },

        {"ares_HI_2.0_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_LE_platinum_15000.h5",{false,0},	15000/5.,1.27,2.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-1.0_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_platinum_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-5.0_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI3_LE_platinum_14995.h5",{false,0},	14995/5.,1.27,-5.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_0.5_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_LE_platinum_15000.h5",{false,0},15000/5.,1.27,0.5,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-3.0_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_LE_platinum_14990.h5",{false,0},	14990/5.,1.27,-3.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES HOLEICE SETS
        * HI1 = 2.0
        * HI2 = -1.0
        * HI3 = -5.0
        * HI4 = 0.5
        * HI5 = -3.0
        */ 
        {"ares_HI_2.0_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_HE_golden_99510.h5",{false,0},	99510/5.,1.27,2.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-1.0_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_golden_98765.h5",{false,0},	98765/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-5.0_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI3_HE_golden_97840.h5",{false,0},	97840/5.,1.27,-5.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_0.5_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_HE_golden_95655.h5",{false,0},	95655/5.,1.27,0.5,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-3.0_golden_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI5_HE_golden_95265.h5",{false,0},	95265/5.,1.27,-3.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },

        {"ares_HI_2.0_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_LE_golden_15000.h5",{false,0},	15000/5.,1.27,2.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-1.0_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_golden_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-5.0_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI3_LE_golden_14995.h5",{false,0},	14995/5.,1.27,-5.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_0.5_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_LE_golden_15000.h5",{false,0},	15000/5.,1.27,0.5,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-3.0_golden_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI5_LE_golden_14990.h5",{false,0},	14990/5.,1.27,-3.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        /*
        * OFFICIAL ARES HOLEICE SETS
        * HI1 = 2.0
        * HI2 = -1.0
        * HI3 = -5.0
        * HI4 = 0.5
        * HI5 = -3.0
        */ 
        {"ares_HI_2.0_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_HE_diamond_99510.h5",{false,0},	99510/5.,1.27,2.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-1.0_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_HE_diamond_99630.h5",{false,0},	99630/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-5.0_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI3_HE_diamond_97840.h5",{false,0},	97840/5.,1.27,-5.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_0.5_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_HE_diamond_95655.h5",{false,0},	95655/5.,1.27,0.5,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_HI_-3.0_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI5_HE_diamond_95265.h5",{false,0},	95265/5.,1.27,-3.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },

        {"ares_HI_2.0_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_LE_diamond_15000.h5",{false,0},	15000/5.,1.27,2.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-1.0_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_Nominal_LE_diamond_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-5.0_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI3_LE_diamond_14995.h5",{false,0},	14995/5.,1.27,-5.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_0.5_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI4_LE_diamond_15000.h5",{false,0},	15000/5.,1.27,0.5,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_HI_-3.0_diamond_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI5_LE_diamond_14990.h5",{false,0},	14990/5.,1.27,-3.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },

        /*
        * OFFICIAL ARES BulkIce Sets
        * BI1 = +5% scattering
        * BI2 = +5% absorption
        * BI4 = -2.5% scatterig -2.5% absorption
        * BI5 = -5% scatterig +5% absorption
        * BI6 = -5% scatterig -5% absorption
        */ 
        {"ares_BI1_diamond_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_HI1_HE_diamond_99510.h5",{false,0},	99510/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI1_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI1_HE_platinum_99810.h5",{false,0},	99810/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI2_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI2_HE_platinum_90260.h5",{false,0},	90260/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI2_platinum_feather_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"BI2_Nominal_HE_platinum_10000.h5",{false,0},	10000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI2_platinum_feather_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI2_LE_platinum_14990.h5",{false,0},	14990/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_BI4_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI4_HE_platinum_30220.h5",{false,0},	30220/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI5_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI5_HE_platinum_37525.h5",{false,0},	37525/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI6_platinum_HE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI6_HE_platinum_36490.h5",{false,0},	36490/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_Generation_data.lic",true)))
        },
        {"ares_BI1_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI1_LE_platinum_14990.h5",{false,0},	14990/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_BI2_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI2_LE_platinum_14990.h5",{false,0},	14990/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_BI4_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI4_LE_platinum_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_BI5_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI5_LE_platinum_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
        {"ares_BI6_platinum_LE",MCSet(MCType::LeptonInjector,mc_dataPath,"Ares_BI6_LE_platinum_15000.h5",{false,0},	15000/5.,1.27,-1.0,1.0,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Ares_LE_Generation_data.lic",true)))
        },
      };
      return simInfo;
    }

    std::vector<MCSet> GetSimulationSets(std::string simulation_tag, std::string mc_dataPath){
      auto simulations = GetSimInfo(mc_dataPath);

      // user specifies directly the simulation they want
      if(not (simulations.find(simulation_tag) == simulations.end()))
        return {simulations.at(simulation_tag)};

      // else we use LE and HE sets
      if((not (simulations.find(simulation_tag+"_LE") == simulations.end())) and (not (simulations.find(simulation_tag+"_LE") == simulations.end())))
        return {simulations.at(simulation_tag+"_LE"),simulations.at(simulation_tag+"_HE")};
      else if ((not (simulations.find(simulation_tag+"_LE") == simulations.end())))
        return {simulations.at(simulation_tag+"_LE")};
      else if ((not (simulations.find(simulation_tag+"_HE") == simulations.end())))
        return {simulations.at(simulation_tag+"_HE")};
      else
        throw std::runtime_error("GolemFit::Invalid simulation tag: "+ simulation_tag + ". Do not append HE or LE postfix to select them both.");
    }

  } // close sterile namespace
/*
  namespace hese {
    std::map<std::string,MCSet> GetSimInfo(std::string mc_dataPath)
    {
       std::map<std::string,MCSet> simInfo = {
        {"nue_scattering",MCSet(MCType::NuGen,mc_dataPath,"nue_scattering.h5",8000*10000.,0.99,0.,1.0,{})},
        {"numu_scattering",MCSet(MCType::NuGen,mc_dataPath,"numu_scattering.h5",2.*8000*10000.,0.99,0.,1.0,{})},
        {"nutau_scattering",MCSet(MCType::NuGen,mc_dataPath,"nutau_scattering.h5",8000*10000.,0.99,0.,1.0,{})},
        {"nue_absorption",MCSet(MCType::NuGen,mc_dataPath,"nue_absorption.h5",8000*10000.,0.99,0.,1.0,{})},
        {"numu_absorption",MCSet(MCType::NuGen,mc_dataPath,"numu_absorption.h5",2.*8000*10000.,0.99,0.,1.0,{})},
        {"nutau_absorption",MCSet(MCType::NuGen,mc_dataPath,"nutau_absorption.h5",8000*10000.,0.99,0.,1.0,{})},
        {"nue_scattering_absorption",MCSet(MCType::NuGen,mc_dataPath,"nue_scattering_absorption.h5",8000*10000.,0.99,0.,1.0,{})},
        {"numu_scattering_absorption",MCSet(MCType::NuGen,mc_dataPath,"numu_scattering_absorption.h5",2.*8000*10000.,0.99,0.,1.0,{})},
        {"nutau_scattering_absorption",MCSet(MCType::NuGen,mc_dataPath,"nutau_scattering_absorption.h5",8000*10000.,0.99,0.,1.0,{})},
        {"nue_dom_081",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_081.h5",8000*10000.,0.81,0.,1.0,{})},
        {"numu_dom_081",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_081.h5",2.*8000*10000.,0.81,0.,1.0,{})},
        {"nutau_dom_081",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_081.h5",8000*10000.,0.81,0.,1.0,{})},
        {"nue_dom_090",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_090.h5",8000*10000.,0.90,0.,1.0,{})},
        {"numu_dom_090",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_090.h5",2.*8000*10000.,0.90,0.,1.0,{})},
        {"nutau_dom_090",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_090.h5",8000*10000.,0.90,0.,1.0,{})},
        {"nue_dom_095",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_095.h5",8000*10000.,0.95,0.,1.0,{})},
        {"numu_dom_095",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_095.h5",2.*8000*10000.,0.95,0.,1.0,{})},
        {"nutau_dom_095",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_095.h5",8000*10000.,0.95,0.,1.0,{})},
        {"nue_dom_108",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_108.h5",8000*10000.,1.08,0.,1.0,{})},
        {"numu_dom_108",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_108.h5",2.*8000*10000.,1.08,0.,1.0,{})},
        {"nutau_dom_108",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_108.h5",8000*10000.,1.08,0.,1.0,{})},
        {"nue_dom_117",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_117.h5",8000*10000.,1.17,0.,1.0,{})},
        {"numu_dom_117",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_117.h5",2.*8000*10000.,1.17,0.,1.0,{})},
        {"nutau_dom_117",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_117.h5",8000*10000.,1.17,0.,1.0,{})},
        {"nue_p2_-3",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_-3.h5",8000*10000.,0.99,-3.,1.0,{})},
        {"numu_p2_-3",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_-3.h5",2.*8000*10000.,0.99,-3.,1.0,{})},
        {"nutau_p2_-3",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_-3.h5",8000*10000.,0.99,-3.,1.0,{})},
        {"nue_p2_-1",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_-1.h5",8000*10000.,0.99,-1.,1.0,{})},
        {"numu_p2_-1",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_-1.h5",2.*8000*10000.,0.99,-1.,1.0,{})},
        {"nutau_p2_-1",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_-1.h5",8000*10000.,0.99,-1.,1.0,{})},
        {"nue_p2_0",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_0.h5",8000*10000.,0.99,0.0,1.0,{})},
        {"numu_p2_0",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_0.h5",2*8000*10000.,0.99,0.0,1.0,{})},
        {"nutau_p2_0",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_0.h5",8000*10000.,0.99,0.0,1.0,{})},
        {"nue_p2_0_half1",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_0_half1.h5",4000*10000.,0.99,0.0,1.0,{})},
        {"numu_p2_0_half1",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_0_half1.h5",2*4000*10000.,0.99,0.0,1.0,{})},
        {"nutau_p2_0_half1",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_0_half1.h5",4000*10000.,0.99,0.0,1.0,{})},
        {"nue_p2_0_half2",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_0_half2.h5",4000*10000.,0.99,0.0,1.0,{})},
        {"numu_p2_0_half2",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_0_half2.h5",2*4000*10000.,0.99,0.0,1.0,{})},
        {"nutau_p2_0_half2",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_0_half2.h5",4000*10000.,0.99,0.0,1.0,{})},
        {"nue_p2_1",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_1.h5",8000*10000.,0.99,1.0,1.0,{})},
        {"numu_p2_1",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_1.h5",2.*8000*10000.,0.99,1.0,1.0,{})},
        {"nutau_p2_1",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_1.h5",8000*10000.,0.99,1.0,1.0,{})},
        {"muongun",MCSet(MCType::MuonGun,mc_dataPath,"muongun.h5",1.0,1.0,0.0,1.0,{})},
      };
      return simInfo;
    }

    std::vector<MCSet> GetSimulationSets(std::string simulation_tag, std::string mc_dataPath){
      auto simulations = GetSimInfo(mc_dataPath);
      if(!(simulation_tag == "absorption" or simulation_tag == "scattering_absorption" or simulation_tag == "scattering" or
         simulation_tag == "dom_081" or simulation_tag == "dom_090" or simulation_tag == "dom_095" or simulation_tag == "dom_108" or simulation_tag == "dom_117" or
         simulation_tag == "p2_-3" or simulation_tag == "p2_-1" or simulation_tag == "p2_0" or simulation_tag == "p2_1" or simulation_tag == "p2_0_half1" or simulation_tag == "p2_0_half2"))
        throw std::runtime_error("GolemFit::GetSimulationSets: simulation_tag :" + simulation_tag + " not registered");
      return {(simulations.at("nue_"+simulation_tag)),
              (simulations.at("numu_"+simulation_tag)),
              (simulations.at("nutau_"+simulation_tag)),
              (simulations.at("muongun"))};
    }
  } // close hese namespace

  namespace magictau {
    std::map<std::string,MCSet> GetSimInfo(std::string mc_dataPath)
    {
       std::map<std::string,MCSet> simInfo = {
        {"nue_baseline",MCSet(MCType::NuGen,mc_dataPath,"nue_baseline.h5",10000*50000.,1.,-1,1.0,{})},
        {"numu_baseline",MCSet(MCType::NuGen,mc_dataPath,"numu_baseline.h5",10000*50000.,1.,-1,1.0,{})},
        {"nutau_baseline",MCSet(MCType::NuGen,mc_dataPath,"nutau_baseline.h5",10000*50000.,1.,-1,1.0,{})},
        {"nue_new",MCSet(MCType::NuGen,mc_dataPath,"nue_new.h5",994*10000.,0.99,0,1.0,{})},
        {"numu_new",MCSet(MCType::NuGen,mc_dataPath,"numu_new.h5",996*10000.*2.,0.99,0,1.0,{})},
        {"nutau_new",MCSet(MCType::NuGen,mc_dataPath,"nutau_new.h5",1000*10000.,0.99,0,1.0,{})},
        {"nue_scattering",MCSet(MCType::NuGen,mc_dataPath,"nue_scattering.h5",8000*10000.,1.0,-1,1.0,{})},
        {"numu_scattering",MCSet(MCType::NuGen,mc_dataPath,"numu_scattering.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nutau_scattering",MCSet(MCType::NuGen,mc_dataPath,"nutau_scattering.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nue_absorption",MCSet(MCType::NuGen,mc_dataPath,"nue_absorption.h5",8000*10000.,1.0,-1,1.0,{})},
        {"numu_absorption",MCSet(MCType::NuGen,mc_dataPath,"numu_absorption.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nutau_absorption",MCSet(MCType::NuGen,mc_dataPath,"nutau_absorption.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nue_scattering_absorption",MCSet(MCType::NuGen,mc_dataPath,"nue_scattering_absorption.h5",8000*10000.,1.0,-1,1.0,{})},
        {"numu_scattering_absorption",MCSet(MCType::NuGen,mc_dataPath,"numu_scattering_absorption.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nutau_scattering_absorption",MCSet(MCType::NuGen,mc_dataPath,"nutau_scattering_absorption.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nue_dom_081",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_081.h5",8000*10000.,0.81,-1,1.0,{})},
        {"numu_dom_081",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_081.h5",8000*10000.,0.81,-1,1.0,{})},
        {"nutau_dom_081",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_081.h5",8000*10000.,0.81,-1,1.0,{})},
        {"nue_dom_090",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_090.h5",8000*10000.,0.90,-1,1.0,{})},
        {"numu_dom_090",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_090.h5",8000*10000.,0.90,-1,1.0,{})},
        {"nutau_dom_090",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_090.h5",8000*10000.,0.90,-1,1.0,{})},
        {"nue_dom_095",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_095.h5",8000*10000.,0.95,-1,1.0,{})},
        {"numu_dom_095",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_095.h5",8000*10000.,0.95,-1,1.0,{})},
        {"nutau_dom_095",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_095.h5",8000*10000.,0.95,-1,1.0,{})},
        {"nue_dom_108",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_108.h5",8000*10000.,1.08,-1,1.0,{})},
        {"numu_dom_108",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_108.h5",8000*10000.,1.08,-1,1.0,{})},
        {"nutau_dom_108",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_108.h5",8000*10000.,1.08,-1,1.0,{})},
        {"nue_dom_117",MCSet(MCType::NuGen,mc_dataPath,"nue_dom_117.h5",8000*10000.,1.17,-1,1.0,{})},
        {"numu_dom_117",MCSet(MCType::NuGen,mc_dataPath,"numu_dom_117.h5",8000*10000.,1.17,-1,1.0,{})},
        {"nutau_dom_117",MCSet(MCType::NuGen,mc_dataPath,"nutau_dom_117.h5",8000*10000.,1.17,-1,1.0,{})},
        {"nue_p2_-3",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_-3.h5",8000*10000.,1.0,-3,1.0,{})},
        {"numu_p2_-3",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_-3.h5",8000*10000.,1.0,-3,1.0,{})},
        {"nutau_p2_-3",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_-3.h5",8000*10000.,1.0,-3,1.0,{})},
        {"nue_p2_-1",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_-1.h5",8000*10000.,1.0,-1,1.0,{})},
        {"numu_p2_-1",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_-1.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nutau_p2_-1",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_-1.h5",8000*10000.,1.0,-1,1.0,{})},
        {"nue_p2_0",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_0.h5",8000*10000.,0.99,0.0,1.0,{})},
        {"numu_p2_0",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_0.h5",2*8000*10000.,0.99,0.0,1.0,{})},
        {"nutau_p2_0",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_0.h5",8000*10000.,0.99,0.0,1.0,{})},
        {"nue_p2_1",MCSet(MCType::NuGen,mc_dataPath,"nue_p2_-1.h5",8000*10000.,1.0,1.0,1.0,{})},
        {"numu_p2_1",MCSet(MCType::NuGen,mc_dataPath,"numu_p2_-1.h5",8000*10000.,1.0,1.0,1.0,{})},
        {"nutau_p2_1",MCSet(MCType::NuGen,mc_dataPath,"nutau_p2_-1.h5",8000*10000.,1.0,1.0,1.0,{})},
        {"muongun_me_tau",MCSet(MCType::MuonGun,mc_dataPath,"muongun_me_tau.h5",1.0,1.0,0.0,1.0,{})},
        {"muongun_he_tau",MCSet(MCType::MuonGun,mc_dataPath,"muongun_he_tau.h5",1.0,1.0,0.0,1.0,{})},
      };
      return simInfo;
    }

    std::vector<MCSet> GetSimulationSets(std::string simulation_tag, std::string mc_dataPath){
      auto simulations = GetSimInfo(mc_dataPath);
      if(!(simulation_tag == "baseline" or simulation_tag == "new" or
         simulation_tag == "absorption" or simulation_tag == "scattering_absorption" or simulation_tag == "scattering" or
         simulation_tag == "dom_081" or simulation_tag == "dom_090" or simulation_tag == "dom_095" or simulation_tag == "dom_108" or simulation_tag == "dom_117" or
         simulation_tag == "p2_-3" or simulation_tag == "p2_-1" or simulation_tag == "p2_0" or simulation_tag == "p2_1"))
        throw std::runtime_error("GolemFit::GetSimulationSets: simulation_tag :" + simulation_tag + " not registered");
      return {(simulations.at("nue_"+simulation_tag)),
              (simulations.at("numu_"+simulation_tag)),
              (simulations.at("nutau_"+simulation_tag)),
              (simulations.at("muongun_me_tau")),
              (simulations.at("muongun_he_tau"))};
    }
  } // close hese namespace
*/

} // close golemfit namespace

#endif
