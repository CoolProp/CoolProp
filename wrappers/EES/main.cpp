//============================================================================================//
//                                                                                            //
//                                  EES - CoolProp interface                                  //
//                                  -------------------------                                 //
//                                                                                            //
//  This dll is an interface between EES and CoolProp.                                        //
//  In EES, external functions need to be implemented in dynamic librairies.  The first       //
//  argument sent by EES is a 256 characters char variable. The second argument is pointer    //
//  structure conaining "double" values.  The third argument is a linked list for the         //
//  input data                                                                                //
//                                                                                            //
//  The arguments are defined as follows :                                                    //
//  - The string variable contains the the definition of the fluids and of their              //
//    concentrations with the input strings concatenated to the fluid name joined by |        //
//    (e.g. "R134a|T|P|D" or "REFPROP-R134a|O|T|P" or                                         //
//    "REFPROP-MIX:R32[0.697615]&R125[0.302385]|V|P|H" (R410A))                               //
//  - mode, which is -1 if to return a default form of the call as string, normal mode        //
//    otherwise                                                                               //
//  - The last value is a linked list of the input values                                     //
//																							  //
//  The file needs to be built in coolprop_ees.dlf, which is the standard extension           //
//  for EES external functions.  If CoolProp has been built to the static library             //
//  CoolPropStaticLibrary.lib, you can build (with visual studio) CoolProp_EES.dlf with:      //
//                                                                                            //
//     link /DEBUG /DLL main.obj CoolPropStaticLibrary.lib /OUT:COOLPROP_EES.dlf              //
//																							  //
//  Only one unit system is used (modified SI - see help). Future versions might              //
//  include a detection of EES current unit system and its definition in the dll              //
//																							  //
//  Ian Bell                                                                                  //
//  Thermodynamics Laboratory                                                                 //
//  University of Liège                                                                       //
//                                                                                            //
//  January 2013                                                                              //
//============================================================================================//

#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <string>
#include <stdio.h> 
#include <string.h>
#include <vector>
#include "CoolProp.h"
#include "CoolPropTools.h"

static const bool EES_DEBUG = false;

// Structure for handling ees calling syntax
struct EesParamRec {
  double value;
  struct EesParamRec *next;
};

// Tell C++ to use the "C" style calling conventions rather than the C++ mangled names

extern "C"  
{
	__declspec (dllexport)  double COOLPROP_EES(char fluid[256], int mode, struct EesParamRec *input_rec)
	{        
		double In1 = _HUGE, In2 = _HUGE, out;  // Two inputs, one output
		int NInputs;           // Ninputs is the number of inputs
		char NInputs_string[3], err_str[1000];
		std::string fluid_string = fluid;
        
		std::string ErrorMsg, Outstr, In1str, In2str, Fluidstr;
		std::vector<std::string> fluid_split;

		if (mode==-1) {	
			strcpy(fluid,"T = coolprop('D','P',101.325,'Q',0,'R134a')"); 
			return 0;
		}

		// Split the string that is passed in at the '~' delimiter that was used to join it
		fluid_split = strsplit(fluid_string,'~');
		if (fluid_split.size() != 4) 
		{
			sprintf(err_str,"fluid[%s] length[%d] not 4 elements long",fluid_string.c_str(),fluid_split.size()); 
			strcpy(fluid,err_str);
            if (EES_DEBUG)
            {
                FILE *fp;
                fp = fopen("log.txt","a+");
                fprintf(fp,"%s %s %g %s %g %s\n",Outstr.c_str(),In1str.c_str(),In1,In2str.c_str(),In2,Fluidstr.c_str());
                fprintf(fp,"%s\n",err_str);
                fclose(fp);
            }
			return 0;
		}
		else{
			Fluidstr = upper(fluid_split[0]);
			Outstr = upper(fluid_split[1]);
			In1str = upper(fluid_split[2]);
			In2str = upper(fluid_split[3]);
		}
		
		// Check the number of inputs
		NInputs = 0;
		EesParamRec * aninput_rec = input_rec;
		while (aninput_rec != 0)
		{
			aninput_rec = aninput_rec->next;
			NInputs++;
		};
		if (NInputs != 2) {
			sprintf(NInputs_string,"%d", NInputs);
			strcpy(fluid,NInputs_string);
			return 0;
		}

		// Get the inputs from the pointer structure sent by EES:
		In1= input_rec->value;
		input_rec=input_rec->next;
		In2=input_rec->value;
		
		//This block can be used to debug the code by writing output or intermediate values to a text file

		if (EES_DEBUG)
        {
            FILE *fp;
            fp = fopen("log.txt","a+");
            fprintf(fp,"%s %s %g %s %g %s\n",Outstr.c_str(),In1str.c_str(),In1,In2str.c_str(),In2,Fluidstr.c_str());
            fclose(fp);
        }

		if (EES_DEBUG)
        {
            // This redirect standard output to file2.txt
            freopen("log_stdout.txt", "w", stdout);
            set_debug_level(10); // Maximum debugging
        }

		try
		{
			out = PropsS(Outstr.c_str(), In1str.c_str(), In1, In2str.c_str(), In2, Fluidstr.c_str());
		}
		catch(...)
		{
			std::string err_str = format("Uncaught error: \"%s\",\"%s\",%g,\"%s\",%g,\"%s\"\n",Outstr.c_str(),In1str.c_str(),In1,In2str.c_str(),In2,Fluidstr.c_str());
			strcpy(fluid, err_str.c_str());
			return 0.0;
		}

		if (!ValidNumber(out))
		{
            std::string err_str = CoolProp::get_global_param_string("errstring");
            // There was an error
            if (EES_DEBUG)
            {
                FILE *fp;
                fp = fopen("log.txt","a+");
                fprintf(fp,"Error: %s \n",err_str.c_str());
                fclose(fp);
            }
			strcpy(fluid,err_str.c_str());
			return 0.0;
		}
		else
        {
            // Check if there was a warning
            std::string warn_string = CoolProp::get_global_param_string("warnstring");
            if (!warn_string.empty())
            {
                if (EES_DEBUG)
                {
                    FILE *fp;
                    fp = fopen("log.txt","a+");
                    fprintf(fp,"Warning: %s \n",warn_string.c_str());
                    fclose(fp);
                }
                // There was a warning, write it back
                strcpy(fluid, warn_string.c_str());
            }
			return out;
        }
	}

};


