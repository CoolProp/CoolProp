

import glob, os

#  Build a map from REFPROP name to CAS code
RP2CAS = {}
for file in glob.glob('C:\\Program Files (x86)\\REFPROP\\fluids\\*.fld'):
    lines = open(file, 'r').readlines()
    root, RPFluid = os.path.split(file)

    for line in lines:
        if line.find('CAS number') > -1:
            CAS_number = line.split('!')[0].strip()
            if not CAS_number:
                raise ValueError(file + line)
            RP2CAS[RPFluid.split('.')[0]] = CAS_number
            break

#  Handle pseudo-pure fluids
for file in glob.glob('C:\\Program Files (x86)\\REFPROP\\fluids\\*.ppf'):
    root, RPFluid = os.path.split(file)
    RP2CAS[RPFluid.split('.')[0]] = RPFluid

fluid_lookup = """1BUTENE	butene	1-Butene	419.29		*	*	*	-	-	-	2.59	0.5	1	0.983	1.079	-	-	-	-	-	-	-	-	-
ACETONE	propanone	Acetone	508.1		-	0.5	-	-	-	-	1.46	0.1	0.2	0.16	9.40E-02	-	-	-	-	-	-	1.39E+04	-	-
AIR	air	Air	132.5306		-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
AMMONIA	ammonia	Ammonia	405.4	R-717	X	X	X	-	-	-	N/A	-	-	-	1.0E-01	-	-	-	-	-	-	1.0E+06	3.50E-01	1.60E+00
ARGON	argon	Argon	150.687	R-740	0	0	0	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
BENZENE	benzene	Benzene	562.02		*	*	*	-	-	-	3.65	0.4	0.2	0.318	2.2E-01	1.9E+03	8.4E-05	2.8E-03	6.4E-05	1.3E-03	1.6E-05	-	-	-
BUTANE	butane	n-Butane	425.125	R-600	N/A	3	N/A	-	-	-	2.15	0.5	0.4	0.485	3.52E-01	-	-	-	-	-	-	-	-	-
C12	dodecane	n-Dodecane	658.1		X	X	X	-	-	-	2.19	0.3	0.4	0.452	3.57E-01	-	-	-	-	-	-	-	-	-
C1CC6	methylcyclohexane		CoolProp error: Your fluid name [] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		-	0.639	-	-	-	-	N/A	0.5	0.6	0.392	1.87	-	-	-	-	-	-	-	-	-
C2BUTENE	cis-2-butene	cis-2-Butene	435.75		*	*	*	-	-	-	2.57	0.4	1		1.15	-	-	-	-	-	-	-	-	-
C3CC6	propylcyclohexane		CoolProp error: Your fluid name [] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		*	*	*	-	-	-	-	-	-	-	2.57	-	-	-	-	-	-	-	-	-
C4F10	perfluorobutane		CoolProp error: Your fluid name [] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		6330	8860	12500	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
C5F12	perfluoropentane		CoolProp error: Your fluid name [] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		6510	9160	13300	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
CF3I	trifluoroiodomethane		CoolProp error: Your fluid name [] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		1*	0.4*	0.1*	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
CO	carbon monoxide	CarbonMonoxide	132.86		-	1.6*	-	-	-	-	0.331	0.04	0.03	0.032	2.70E-02	-	-	-	-	-	-	-	-	-
CO2	carbon dioxide	CarbonDioxide	304.1282		1	1	1	-	-	-	0.108	-	-	-	-	-	-	-	-	-	-	-	-	-
COS	carbonyl sulphide	CarbonylSulfide	378.77		97	27	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
CYCLOHEX	cyclohexane	CycloHexane	553.64		X	X	X	-	-	-	N/A	N/A	N/A	N/A	N/A	-	-	-	-	-	-	-	-	-
CYCLOPEN	cyclopentane	Cyclopentane	511.72		*	*	*	-	-	-	N/A	N/A	N/A	N/A	N/A	-	-	-	-	-	-	-	-	-
CYCLOPRO	cyclopropane	CycloPropane	398.3		*	*	*	-	-	-	N/A	N/A	N/A	N/A	N/A	-	-	-	-	-	-	-	-	-
D2	deuterium	Deuterium	CoolProp error: Your fluid name [Deuterium] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		0	0	0	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
D2O	deuterium oxide	DeuteriumOxide	CoolProp error: Your fluid name [DeuteriumOxide] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
D4	octamethylcyclotetrasiloxane	D4	586.5		N/A	N/A	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
D5	decamethylcyclotetrasiloxane	D5	619.15		N/A	N/A	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
D6	dodecamethylcyclotetrasiloxane	D6	645.78		N/A	N/A	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
DECANE	decane	n-Decane	617.7		*	*	*	-	-	-	2.45	0.4	0.5	0.509	3.84E-01	-	-	-	-	-	-	-	-	-
DMC	dimethyl carbonate	DimethylCarbonate	557		N/A	N/A	N/A	-	-	-	N/A	-	-	-	2.50E-02	-	-	-	-	-	-	-	-	-
DME	dimethylether	DimethylEther	400.378		1	1	<<1	-	-	-	1.66	0.3	0.3		1.89E-01	-	-	-	-	-	-	-	-	-
ETHANE	ethane	Ethane	305.322	R-170	N/A	2.9	N/A	-	-	-	1.46	0.1	0.1	0.121	1.23E-01	-	-	-	-	-	-	-	-	-
ETHANOL	ethanol	Ethanol	514.71		-	-	-	-	-	-	1.95	0.2	0.3	0.317	3.99E-01	-	-	-	-	-	-	1.56E+06	-	-
ETHYLENE	ethene	Ethylene	282.35	R-1150	N/A	6.8	N/A	-	-	-	3.45	1	1	1	1.0E+00	6.4E-01	1.4E-11	7.9E-11	9.0E-12	7.1E-11	1.3E-12	-	-	-
FLUORINE	fluorine	Fluorine	CoolProp error: Your fluid name [Fluorine] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		-	-	-	-	-	-	4.86	-	-	-	-	-	-	-	-	-	-	-	-	-
H2S	hydrogen sulfide	HydrogenSulfide	373.1		-	-	-	-	-	-	6.89	-	-	-	-	2.2E-01	-	-	-	-	-	2.3E+09	-	-
HELIUM	helium	Helium	5.1953	R-704	-	-	-	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
HEPTANE	heptane	n-Heptane	540.13		*	*	*	-	-	-	2.58	0.5	0.5	0.592		-	-	-	-	-	-	-	-	-
HEXANE	hexane	n-Hexane	507.82		*	3.1	*	-	-	-	2.57	0.5	0.4	0.495	4.94E-01	-	-	-	-	-	-	-	-	-
HYDROGEN		Hydrogen	33.145	R-702	-	-	-	-	-	-	N/A	-	-	-	4.94E-01	-	-	-	-	-	-	-	-	-
IBUTENE	2-methyl-1-propene/methylpropene/isobutene/isobutylene	Isobutene	418.09		*	*	*	-	-	-	N/A	0.6	0.6		6.27E-01	-	-	-	-	-	-	6.67E+04
IHEXANE	2-methylpentane (methylpentane)	Isohexane	497.7		*	*	*	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
IPENTANE	2-methylbutane	Isopentane	460.35		*	*	*	-	-	-	1.8	0.3	0.3		4.05E-01	-	-	-	-	-	-	-	-	-
ISOBUTAN	2-methylpropane	IsoButane	407.817		*	*	*	-	-	-	1.74	0.4	0.3		3.07E-01	-	-	-	-	-	-	-	-	-
KRYPTON	kr	Krypton	209.48	R-784	-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
MD2M	decamethyltetrasiloxane	MD2M	599.4		*	*	*	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
MD3M	dodecamethylcyclotetrasiloxane	MD3M	628.36		*	*	*	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
MD4M	tetradecamethylhexasiloxane	MD4M	653.2		*	*	*	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
MDM	octamethyltrisiloxane	MDM	564.09		*	*	*	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
METHANE	methane	Methane	190.564		72	25	7.6	-	-	-	2.72	0.007	0.007		6.00E-03	-	-	-	-	-	-	-	-	-
METHANOL	methanol	Methanol	CoolProp error: Your fluid name [Methanol] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		N/A	2.8	N/A	-	-	-	1.44	0.2	0.1	0.178	1.40E-01	-	-	-	-	-	-	1.37E+04	-	-
MLINOLEA	methyl linoleate (methyl (z,z)-9,12-octadecadienoate)	MethylLinoleate	799		N/A	N/A	N/A	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
MLINOLEN	methyl linolenate (methyl (z,z,z)-9,12,15-octadecatrienoate)	MethylLinolenate	772		N/A	N/A	N/A	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
MM	hexamethyldisiloxane	MM	518.75		N/A	N/A	N/A	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
MOLEATE	methyl oleate (methyl cis-9-octadecenoate)	MethylOleate	782		N/A	N/A	N/A	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
MPALMITA	methyl hexadecanoate 	MethylPalmitate	755		N/A	N/A	N/A	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
MSTEARAT	methyl octadecanoate	MethylStearate	775		N/A	N/A	N/A	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
N2O	nitrous oxide	NitrousOxide	309.52		290	320	180	-	-	-	N/A	N/A	N/A	N/A	N/A	-	-	-	-	-	-	-	-	-
NEON	neon	Neon	44.4918	R-720	-	-	-	-	-	-	N/A	N/A	N/A	N/A	N/A	-	-	-	-	-	-	-	-	-
NEOPENTN	neopentane (2,2-dimethylpropane)	Neopentane	433.74		*	*	*	-	-	-	2.25	-	-	-	1.73E-01	-	-	-	-	-	-	-	-	-
NF3	nitrogen trifluoride	NitrogenTrifluoride	CoolProp error: Your fluid name [NitrogenTrifluoride] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		12300	17200	20700	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
NITROGEN	nitrogen	Nitrogen	126.192	R-728	-	-	-	-	-	-	N/A	-	-	-		-	-	-	-	-	-	-	-	-
NONANE	nonane	n-Nonane	594.55		*	*	*	-	-	-	2.29	0.4	0.5	0.463	4.14E-01	-	-	-	-	-	-	-	-	-
OCTANE	octane	n-Octane	569.32		*	*	*	-	-	-	2.41	0.5	0.5	0.544	4.53E-01	-	-	-	-	-	-	-	-	-
ORTHOHYD	orthohydrogen	OrthoHydrogen	33.22		-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
OXYGEN	oxygen	Oxygen	154.581		-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
PARAHYD	parahydrogen	ParaHydrogen	32.938		-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
PENTANE	pentane	n-Pentane	469.7	R-601	*	*	*	-	-	-	N/A	0.3	0.4	0.387	3.95E-01	-	-	-	-	-	-	-	-	-
PROPANE	propane	n-Propane	369.89	R-290	*	3	*	-	-	-	2.24	0.5	0.4	0.518	1.76E-01	-	-	-	-	-	-	-	-	-
PROPYLEN	propylene	Propylene	364.211		*	3.1	*	-	-	-	2.64	0.6	1	1.06	1.12E+00	-	-	-	-	-	-	-	-	-
PROPYNE	propyne/methylacetylene	Propyne	CoolProp error: Your fluid name [Propyne] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid		*	*	*	-	-	-	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A
R11	trichlorofluoromethane 	R11	471.06	CFC-11	6730	4750	1620	1	1	1	541	-	-	-	-	-	-	-	-	-	-	-	-	-
R113	1,1,2-trichloro-1,2,2-trifluoroethane	R113	CoolProp error: Your fluid name [R113] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	CFC-113  	6540	6130	2700	0.8	0.8	1	659	-	-	-	-	-	-	-	-	-	-	-	-	-
R114	1,2-dichloro-1,1,2,2-tetrafluoroethane	R114	CoolProp error: Your fluid name [R114] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	CFC-114	8040	10000	8730	1	1	1	1110	-	-	-	-	-	-	-	-	-	-	-	-	-
R115	chloropentafluoroethane  	R115	CoolProp error: Your fluid name [R115] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	CFC-115  	5310	7370	9990	0.6	0.6	0.44	1080	-	-	-	-	-	-	-	-	-	-	-	-	-
R116	hexafluoroethane/perfluoroethane	R116	293.03	FC-116	N/A	9200-12200 (2003)	N/A	-	-	-	1380	-	-	-	-	-	-	-	-	-	-	-	-	-
R12	dichlorodifluoromethane  	R12	CoolProp error: Your fluid name [R12] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	CFC-12	11000	10900	5200	1	1	1	1040	-	-	-	-	-	-	-	-	-	-	-	-	-
R123	2,2-dichloro-1,1,1-trifluoroethane 	R123	456.82	HCFC-123	273	77	24	0.02	0.02	0.02	12.3	-	-	-	-	-	-	-	-	-	-	-	-	-
R1234YF	2,3,3,3-tetrafluoroprop-1-ene	R1234yf	367.85	R-1234yf/HFO-1234yf	N/A	4	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R1234ZE	trans-1,3,3,3-tetrafluoropropene 	R1234ze(E)	382.52	R-1234ze	N/A	6	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R124	1-chloro-1,2,2,2-tetrafluoroethane	R124	395.425	HCFC-124	2070	609	185	0.022	0.022	0.022	55.3	-	-	-	-	-	-	-	-	-	-	-	-	-
R125	pentafluoroethane	R125	339.173	HFC-125  	6350	3500	1100	-	-	-	354	-	-	-	-	-	-	-	-	-	-	-	-	-
R13	chlorotrifluoromethane	R13	CoolProp error: Your fluid name [R13] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	CFC-13 	10800	14400	16400	1	1	-	1390	-	-	-	-	-	-	-	-	-	-	-	-	-
R134A	1,1,1,2-tetrafluoroethane	R134a	374.21	HFC-134a	3830	1430	435	-	-	-	144	-	-	-	-	-	-	-	-	-	-	-	-	-
R14	tetrafluoromethane/perfluoromethane	R14	CoolProp error: Your fluid name [R14] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	FC-14	5210	7390	11200	-	-	-	697	-	-	-	-	-	-	-	-	-	-	-	-	-
R141B	1,1-dichloro-1-fluoroethane 	R141b	CoolProp error: Your fluid name [R141B] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	HCFC-141b	2250	725	220	0.11	0.11	0.12	80.6	-	-	-	-	-	-	-	-	-	-	-	-	-
R142B	1-chloro-1,1-difluoroethane 	R142b	CoolProp error: Your fluid name [R142B] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	HCFC-142b	5490	2310	705	0.065	0.065	0.07	228	-	-	-	-	-	-	-	-	-	-	-	-	-
R143A	1,1,1-trifluoroethane	R143a	345.857	HFC-143a 	5890	4470	1590	-	-	-	487	-	-	-	-	-	-	-	-	-	-	-	-	-
R152A	1,1-difluoroethane	R152A	386.411	HFC-152a 	437	124	38	-	-	-	15.5	-	-	-	-	-	-	-	-	-	-	-	-	-
R161	fluoroethane/ethylfluoride	R161	375.25		N/A	10	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R21	dichlorofluoromethane 	R21	CoolProp error: Your fluid name [R21] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	HCFC-21 	530	151	46	0.04	0.04	N/A	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R218	octafluoropropane/perfluoropropane	R218	345.02	F-218	6310	8830	12500	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R22	chlorodifluoromethane 	R22	369.295	HCFC-22	5160	1810	549	0.055	0.055	0.05	194	-	-	-	-	-	-	-	-	-	-	-	-	-
R227EA	1,1,1,2,3,3,3-heptafluoropropane	R227EA	374.9	HFC-227ea 	5310	3220	1040	-	-	-	365	-	-	-	-	-	-	-	-	-	-	-	-	-
R23	trifluoromethane	R23	299.293	HFC-23 	12000	14800	12200	-	-	-	1340	-	-	-	-	-	-	-	-	-	-	-	-	-
R236EA	1,1,1,2,3,3-hexafluoropropane 	R236EA	412.44	HFC-236ea	N/A	1200	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R236FA	1,1,1,3,3,3-hexafluoropropane	R236FA	398.07	HFC-236fa 	8100	9810	7660	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R245CA	1,1,2,2,3-pentafluoropropane 	R245CA	CoolProp error: Your fluid name [R245CA] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	HFC-245ca 	N/A	N/A	N/A	-	-	-	67.5	-	-	-	-	-	-	-	-	-	-	-	-	-
R245FA	1,1,1,3,3-pentafluoropropane	R245fa	427.16	HFC-245fa 	3380	1030	314	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R32	difluoromethane	R32	351.255	HFC-32	2330	675	205	-	-	-	64.2	-	-	-	-	-	-	-	-	-	-	-	-	-
R365MFC	1,1,1,3,3-pentafluorobutane  	R365MFC	460	HFC-365mfc	2520	794	241	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R404A	44% r-125, 52% r143a, r134a	R404A	345.27	Blend	N/A	3900	N/A	0.04	0	0	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R407C	23% r-32, 25% r-125, 52% r134a	R407C	359.345	Blend	N/A	1800	N/A	0	N/A	N/A	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R41	fluoromethane	R41	317.28	HFC-41 	323	92	28	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R410A	50% r-32, 50% r-125	R410A	344.494	Blend	N/A	2088	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
R507A	50% r-125, 50% r-143a	R507A	343.765	Blend	N/A	3985	N/A	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
RC318	octafluorocyclobutane/perfluorocyclobutane	RC318	CoolProp error: Your fluid name [RC318] is not a CoolProp fluid, a REFPROP fluid, a brine or a liquid	FC-C318 	7310	10300	14700	-	-	-	1010	-	-	-	-	-	-	-	-	-	-	-	-	-
SF6	sulfur hexafluoride/sulphur hexafluoride	SulfurHexafluoride	318.7232		16300	22800	32600	-	-	-	2760	-	-	-	-	-	-	-	-	-	-	-	-	-
SO2	sulfur dioxide/sulphur dioxide	SulfurDioxide	430.64	R-764	**	**	**	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
T2BUTENE	trans-2-butene	trans-2-Butene	428.61	C4H8	*	*	*	-	-	-	2.57	-	-	-	1.13E+00	-	-	-	-	-	-	-	-	-
TOLUENE	methylbenzane	Toluene	591.75		N/A	3.3	N/A	-	-	-	1.95	0.5	0.6	0.565	6.4E-01	3.3E-01	7.0E-05	7.0E-04	5.0E-05	5.8E-04	1.6E-05	2.6E+05	-	-
WATER		Water	647.096		***	***	***	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-
XENON		Xenon	289.733		-	-	-	-	-	-	N/A	-	-	-	-	-	-	-	-	-	-	-	-	-"""

name_dict, ODP_dict, GWP20_dict, GWP100_dict, GWP500_dict = {}, {}, {}, {}, {}

for row in fluid_lookup.split('\n'):
    a = row.split('\t')
    #  Refprop fluid name
    RPName = a[0].strip()

    #  CAS number for this fluid
    CAS = RP2CAS[RPName]

    name_dict[CAS] = a[2].strip()
    ODP_dict[CAS] = a[10].strip()
    GWP20_dict[CAS] = a[5].strip()
    GWP100_dict[CAS] = a[6].strip()
    GWP500_dict[CAS] = a[7].strip()

ASHRAE34data = """R11	A1
R12	A1
R13	A1
R21	B1
R22	A1
R23	A1
R30	B2
R32	A2
R40	B2
METHANE	A3
R113	A1
R114	A1
R115	A1
R116	A1
R123	B1
R124	A1
R125F	A1
R134A	A1
R142B	A2
R143A	A2
R152A	A2
ETHANE	A3
DME	A3
R218	A1
R227EA	A1
R236FA	A1
R245FA	B1
PROPANE	A3
RC318	A1
BUTANE	A3
ISOBUTANE	A3
IPENTANE	A3
HYDROGEN	A3
HELIUM	A1
AMMONIA	B2
WATER	A1
NEON	A1
NITROGEN	A1
ARGON	A1
CO2	A1
SO2	B1
ETHYLENE	A3
PROPYLEN	A3
R404A	A1
R507A	A1
R410A	A1
R407C	A1
R1234YF	A2L"""

ASHRAE34_dict = {}
for row in ASHRAE34data.split('\n'):
    a = row.split('\t')
    if a[0] in RP2CAS:
        ASHRAE34_dict[RP2CAS[a[0]]] = a[1]
    else:
        print('Missing CAS number for ' + a[0])

fluids = """:'1BUTENE.FLD','ACETONE.FLD','AIR.PPF','AMMONIA.FLD','ARGON.FLD',
:'BENZENE.FLD','BUTANE.FLD','C1CC6.FLD','C2BUTENE.FLD','C3CC6.FLD',
:'C4F10.FLD','C5F12.FLD','C12.FLD','CF3I.FLD','CO.FLD','CO2.FLD',
:'COS.FLD','CYCLOHEX.FLD','CYCLOPEN.FLD','CYCLOPRO.FLD','D2.FLD',
:'D2O.FLD','D4.FLD','D5.FLD','D6.FLD','DECANE.FLD','DMC.FLD',
:'DME.FLD','ETHANE.FLD','ETHANOL.FLD','ETHYLENE.FLD','FLUORINE.FLD'
:,'H2S.FLD','HELIUM.FLD','HEPTANE.FLD','HEXANE.FLD','HYDROGEN.FLD',
:'IBUTENE.FLD','IHEXANE.FLD','IPENTANE.FLD','ISOBUTAN.FLD',
:'KRYPTON.FLD','MD2M.FLD','MD3M.FLD','MDM.FLD','METHANE.FLD',
:'METHANOL.FLD','MLINOLEA.FLD','MLINOLEN.FLD','MM.FLD',
:'MOLEATE.FLD','MPALMITA.FLD','MSTEARAT.FLD','N2O.FLD','NEON.FLD',
:'NEOPENTN.FLD','NF3.FLD','NITROGEN.FLD','NONANE.FLD',
:'OCTANE.FLD','ORTHOHYD.FLD','OXYGEN.FLD','PARAHYD.FLD',
:'PENTANE.FLD','PROPANE.FLD','PROPYLEN.FLD','PROPYNE.FLD',
:'R32.FLD','R41.FLD','R115.FLD','R116.FLD','R124.FLD','R125.FLD',
:'R141B.FLD','R142B.FLD','R143A.FLD','R161.FLD','R218.FLD',
:'R227EA.FLD','R236EA.FLD','R236FA.FLD','R245CA.FLD','R245FA.FLD',
:'R365MFC.FLD','R507A.PPF','R1234YF.FLD','R1234ZE.FLD','SF6.FLD',
:'SO2.FLD','T2BUTENE.FLD','TOLUENE.FLD','WATER.FLD','XENON.FLD',
:'R11.FLD','R12.FLD','R13.FLD','R14.FLD','R21.FLD','R22.FLD',
:'R23.FLD','R113.FLD','R114.FLD','R123.FLD','R134A.FLD','R152A.FLD',
:'R404A.PPF','R407C.PPF','R410A.PPF','RC318.FLD'"""

HH = """:'0','2','0','3','0',
:'2','1','2','0','NA',
:'1','NA','2','1','1','1',
:'3','1','2','2','NA',
:'NA','NA','NA','NA','2','2',
:'1','1','2','2','4'
:,'4','0','1','2','0',
:'1','2','1','1',
:'0','1','1','1','0',
:'2','NA','NA','2',
:'2','1','0','1','0',
:'1','1','0','2',
:'2','NA','0','NA',
:'2','1','1','1',
:'1','2','1','1','1','1',
:'1','1','1','NA','2',
:'1','NA','1','NA','2',
:'NA','1','1','1','1',
:'3','0','2','0','0',
:'1','1','NA','NA','NA','1',
:'1','1','1','2','1','1',
:'1','1','1','1'"""

FH = """:'4','3','0','1','0',
:'3','4','3','4','NA',
:'0','NA','2','0','4','0',
:'4','3','3','2','NA',
:'NA','NA','NA','NA','2','3',
:'4','4','3','4','3'
:,'4','0','3','3','4',
:'4','3','4','4',
:'0','2','2','2','4',
:'3','NA','NA','3',
:'1','0','1','0','0',
:'4','0','0','3',
:'3','NA','0','NA',
:'4','4','4','4',
:'4','3','0','0','1','1',
:'1','1','1','NA','1',
:'0','NA','0','NA','0',
:'NA','1','2','2','0',
:'0','4','3','0','0',
:'1','1','NA','NA','NA','1',
:'1','0','0','1','1','4',
:'1','1','1','0' """

PH = """:'0','0','0','0','0',
:'0','0','0','1','NA',
:'0','NA','0','0','3','0',
:'1','0','1','0','NA',
:'NA','NA','NA','NA','0','0',
:'2','0','0','2','0'
:,'0','0','0','0','0',
:'2','1','0','0',
:'0','0','0','0','0',
:'0','NA','NA','1',
:'0','0','0','0','0',
:'0','3','0','0',
:'0','NA','0','NA',
:'0','0','1','1',
:'1','2','1','1','0','0',
:'0','0','0','NA','1',
:'0','NA','1','NA','1',
:'NA','0','0','0','0',
:'0','1','0','0','3',
:'0','0','NA','NA','NA','0',
:'0','0','0','0','0','1',
:'0','0','0','2'"""

pp_fluids = fluids.replace(':', '').replace('\n', '').replace('.FLD', '').replace('.PPF', '').replace("'", '').split(",")
pp_HH = HH.replace(':', '').replace('\n', '').replace("'", '').split(",")
pp_FH = FH.replace(':', '').replace('\n', '').replace("'", '').split(",")
pp_PH = PH.replace(':', '').replace('\n', '').replace("'", '').split(",")

HH_dict = {RP2CAS[k]: v for k, v in zip(pp_fluids, pp_HH)}
FH_dict = {RP2CAS[k]: v for k, v in zip(pp_fluids, pp_FH)}
PH_dict = {RP2CAS[k]: v for k, v in zip(pp_fluids, pp_PH)}


def get_env_data(fluid):
    a = dict(
        HH=HH_dict[fluid],
        FH=FH_dict[fluid],
        PH=PH_dict[fluid],
        ODP=ODP_dict[fluid],
        GWP20=GWP20_dict[fluid],
        GWP100=GWP100_dict[fluid],
        GWP500=GWP500_dict[fluid]
        )

    for k, v in a.iteritems():
        try:
            a[k] = int(v)
        except ValueError:
            try:
                a[k] = float(v)
            except ValueError:
                a[k] = -1

    for term in ['GWP100', 'GWP20', 'GWP500', 'ODP']:
        try:
            a[term] = float(a[term])
        except TypeError:
            a[term] = -1

    if fluid in ASHRAE34_dict:
        a['ASHRAE34'] = ASHRAE34_dict[fluid]
    else:
        a['ASHRAE34'] = 'UNKNOWN'

    a['Name'] = name_dict[fluid]

    return a


import json

code = {}
for fluid in pp_fluids:
    fluid = RP2CAS[fluid]
    if name_dict[fluid]:
        code[fluid] = get_env_data(fluid)
    else:
        continue

f = open('DTU_environmental.json', 'w')
f.write(json.dumps(code, sort_keys=True, indent=2, separators=(',', ': ')))
f.close()
