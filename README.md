# PV-BESS Tool [PVBT] (Analysis and Sizing tool for the small-scale PV/BESS)

![PVBT](https://user-images.githubusercontent.com/69669859/104011567-3a54f780-51a6-11eb-9295-5aebd8e59f5e.png)

[![DOI](https://zenodo.org/badge/327884914.svg)](https://zenodo.org/badge/latestdoi/327884914)

[![View PV-BESS Analysis and Sizing Tool (PVBT) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/98279-pv-bess-analysis-and-sizing-tool-pvbt)

This tool was validated and detailed in the following paper: 

A. A. R. Mohamed, R. J. Best, X. A. Liu and D. J. Morrow, "A Comprehensive Robust Techno-Economic Analysis and Sizing Tool for the Small-Scale PV and BESS," in IEEE Transactions on Energy Conversion, 2021, doi: 10.1109/TEC.2021.3107103.

The full paper link: https://ieeexplore.ieee.org/document/9521731  ||  Open Access version from: https://bit.ly/3hspvsu

The repository contains two programs: 1) PVBT and 2) PVBTOptimization. The first one is for the analysis only, the second one is for the PV/BESS sizing optimization and analysis.

The PVBT tool utilizes a real-time BESS control method that aim to maximize the PV self-consumption and energy arbitrage that has been validated using real measurements in addition to integrating a rigorous ageing model to determine the loss in savings due to the capacity degradation. The PVBT outputs are:
1.	The net household demand with and without the BESS. 
2.	Electricity bill with and without the BESS.
3.	BESS power dispatch.
4.	BESS state of charge.
5.	Battery degradation and SoH. 
6.	PV self-consumption with and without the BESS. 
7.	Self-sufficiency with and without the BESS.
8.	Power curtailed with and without the BESS.
9.	Exported power to the grid with and without the BESS.

In addition to five cost benefit analysis (CBA) are conducted for:

1.	CBA for the PV investment according to the warrantied lifetime.
2.	CBA for the BESS investment according to the warrantied lifetime.
3.	CBA for the BESS investment according to the lifetime based on the minimum state of health (end of lifetime based on SoH).
4.	CBA for the PV and the BESS according to the warrantied lifetime of the PV and the BESS.
5.	CBA for the PV and the BESS according to the PV warrantied lifetime and the BESS lifetime based on the minimum state of health.

The aim of the optimization formulation (PVBTOptimization) is to find the optimal sizes of PV only with or without BESS, BESS only in presence of PV, and PV with BESS sequentially. The optimization objective is to maximize the profitability through maximizing the net present value. 

Please check the PVBTGuide.pdf file for more details and guidance on how to use the code. 

This code has been developed by [Ahmed A.Raouf Mohamed](https://pure.qub.ac.uk/en/persons/ahmed-mohamed) - EPIC Research Cluster, School of Electronics, Electrical Engineering and Computer Science at Queen's University in Belfast, UK. This work is part of [SPIRE 2 Project](https://www.ulster.ac.uk/spire2/the-project). 

For any inquiry: amohamed06@qub.ac.uk / AARaoufM@gmail.com 
[![twitter2](https://user-images.githubusercontent.com/69669859/97111234-a068cd00-16d5-11eb-9559-ff4b8946c0d8.png)](https://twitter.com/RA2OOOF)

v1.0 First release (01/2021).

v1.1 Added "FB" that represents fixed component of the BESS prices which is not affected by the annual declining rate in BESS prices (01/2021).
v1.2 An option to discharge the battery starting from the end of low tariff period (02/2021).
Copyright @ 2021 
